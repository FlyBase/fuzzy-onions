# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
#
# Redistribution and use of this script, with or without modifications,
# is permitted provided that the following conditions are met:
#
# 1. Redistributions of this script must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

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
