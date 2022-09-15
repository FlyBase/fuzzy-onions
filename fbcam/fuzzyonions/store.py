# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import logging

import click
from click_shell.core import make_click_shell

from fbcam.fuzzyonions.scea import FileStore


class Store(object):
    """Represents the local file storage for raw data."""

    def __init__(self, config):
        self._prod = FileStore(config.get('store', 'production'))
        self._staging = FileStore(config.get('store', 'staging'), staging=True)

    @property
    def datasets(self):
        """Gets all the datasets in the store."""

        return self._prod.datasets + self._staging.datasets

    def get(self, dsid, download=True):
        """Gets a single dataset.

        :param dsid: the dataset ID
        :param download: if True, the dataset will be downloaded if it
            is not already available locally
        """

        ds = self._prod.get(dsid, download)
        if ds is None:
            ds = self._staging(dsid, download)
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


@click.group(name="store", invoke_without_command=True)
@click.pass_context
def store(ctx):
    """Manage the local file store."""

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-store> ")
        shell.cmdloop()


@store.command()
@click.argument('dsid')
@click.pass_obj
def download(ctx, dsid):
    """Download a SCEA dataset.

    This command fetches the data for the specified dataset from the
    SCEA server.
    """

    ds = ctx.raw_store.get(dsid)
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
