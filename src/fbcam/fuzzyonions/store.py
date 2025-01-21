# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import logging
from datetime import datetime

import click
from click_shell.core import make_click_shell

from fbcam.fuzzyonions.scea import FileStore


class Store(object):
    """Represents the local file storage for raw data."""

    def __init__(self, config):
        self._prod = FileStore(config.get('store', 'production'))
        self._staging = FileStore(config.get('store', 'staging'), staging=True)
        custom_path = config.get('store', 'custom', fallback=None)
        if custom_path is not None:
            self._custom = FileStore(custom_path)
        else:
            self._custom = None

    @property
    def datasets(self):
        """Gets all the datasets in the store."""

        ds = self._prod.datasets + self._staging.datasets
        if self._custom is not None:
            ds += self._custom.datasets
        return ds

    def get(self, dsid, download=True, staging=False):
        """Gets a single dataset.

        :param dsid: the dataset ID
        :param download: if True, the dataset will be downloaded if it
            is not already available locally
        :param staging: if True, the dataset will be downloaded from
            the staging server, regardless of whether it is available
            on the production server
        """

        # First look up only locally
        ds = self._prod.get(dsid, download=False)
        if ds is None:
            ds = self._staging.get(dsid, download=False)
        if ds is None and self._custom is not None:
            ds = self._custom.get(dsid, download=False)

        # Not available locally, download?
        if ds is None and download:
            if not staging:
                ds = self._prod.get(dsid, download=True)
            if ds is None:
                ds = self._staging.get(dsid, download=True)

        return ds

    def update(self):
        """Update the store from the production server."""

        updated = []
        for dataset in self._staging.datasets:
            if dataset.id not in [ds.id for ds in self._prod.datasets]:
                logging.info(f"Checking for {dataset.id} on production")
                ds = self._prod.get(dataset.id)
                if ds is not None:
                    updated.append(dataset.id)
            else:
                logging.info(f"Removing {dataset.id} from staging")
                updated.append(dataset.id)

        if len(updated) > 0:
            self._staging.delete(updated)

        return updated

    def delete(self, dsids):
        """Delete specified datasets from the store."""

        self._staging.delete(dsids)
        self._prod.delete(dsids)
        if self._custom is not None:
            self._custom.delete(dsids)


@click.group(name="store", invoke_without_command=True)
@click.pass_context
def store(ctx):
    """Manage the local file store."""

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-store> ")
        shell.cmdloop()


@store.command()
@click.argument('dsid')
@click.option(
    '--force',
    '-f',
    is_flag=True,
    default=False,
    help="Download the dataset even if it already exists in the local store.",
)
@click.option(
    '--staging',
    '-s',
    is_flag=True,
    default=False,
    help="Force downloading from the staging server.",
)
@click.pass_obj
def download(ctx, dsid, force, staging):
    """Download a SCEA dataset.

    This command fetches the data for the specified dataset from the
    SCEA server.
    """

    if force:
        ctx.raw_store.delete([dsid])

    ds = ctx.raw_store.get(dsid, staging=staging)
    if ds:
        ctx.tracker.add_dataset(dsid, ds.staging)
        ctx.tracker.save()


@store.command('list')
@click.pass_obj
def list_datasets(ctx):
    """List datasets in the local store."""

    n = len(ctx.raw_store.datasets)
    if n == 0:
        print("No local datasets available.")
        return

    print(f"{n} dataset(s) available:")
    for d in ctx.raw_store.datasets:
        print(f"{d.id}")


@store.command()
@click.pass_obj
def update(ctx):
    """Update all locally available datasets.

    This command checks whether datasets from the SCEA staging server
    have been moved to production and if so, updates the local cache
    with the production files.
    """

    upds = ctx.raw_store.update()
    for dsid in upds:
        ctx.tracker.promote_to_production(dsid)
    if len(upds) > 0:
        ctx.tracker.save()


@store.command()
@click.option(
    '--since',
    '-s',
    type=click.DateTime(),
    default=None,
    help="Ignore datasets older than the specified date.",
)
@click.option(
    '--download',
    '-b',
    is_flag=True,
    default=False,
    help="Download the datasets to the local store.",
)
@click.pass_obj
def findnew(ctx, since, download):
    """Find new datasets on the remote store.

    This command checks whether new datasets have been made available
    on the SCEA staging server.
    """

    experiments = ctx.raw_store._staging.get_experiments_list()
    known_ids = [d.id for d in ctx.raw_store.datasets]
    downloaded = 0

    for experiment in experiments:
        accession = experiment['experimentAccession']
        if accession in known_ids:
            continue

        load_date = datetime.strptime(experiment['loadDate'], '%d-%m-%Y')
        if since is not None and load_date < since:
            continue

        print(f"{accession}\t{load_date:%F}\t{experiment['experimentDescription']}")
        if download:
            ds = ctx.raw_store.get(accession)
            if ds:
                downloaded += 1
                tds = ctx.tracker.add_dataset(accession, ds.staging)

                if (
                    'inferred cell type - authors labels'
                    in experiment['experimentalFactors']
                    or 'inferred cell type - ontology labels'
                    in experiment['experimentalFactors']
                ):
                    tds.cell_types.set_available(date=load_date)

    if downloaded > 0:
        ctx.tracker.save()
