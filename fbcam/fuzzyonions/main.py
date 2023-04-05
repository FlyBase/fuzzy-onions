# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2021 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import logging
from os import getenv
from os.path import exists
from configparser import ConfigParser

import click
from click_shell import shell
from IPython import embed

from fbcam.fuzzyonions import __version__
from fbcam.fuzzyonions.explorer import explorer
from fbcam.fuzzyonions.curation import curate
from fbcam.fuzzyonions.store import store, Store
from fbcam.fuzzyonions.tracker import tracker, DatasetTracker
from fbcam.fuzzyonions.database import DatabaseHelper
from fbcam.fuzzyonions.dsfinder import discover
from fbcam.fuzzyonions.annotations import annots

prog_name = "fuzzyonions"
prog_notice = f"""\
{prog_name} {__version__}
Copyright © 2021 Damien Goutte-Gattat

This program is released under the terms of the 1-clause BSD licence.
"""


class FzoContext(object):
    def __init__(self, config_file):
        self._config_file = config_file
        self._config = ConfigParser()

        self.reset()

    @property
    def has_config(self):
        return self._has_config

    @property
    def config(self):
        return self._config

    @property
    def config_file(self):
        return self._config_file

    def reset(self, options=None):
        self._store = None
        self._tracker = None
        self._database = None

        self._config.clear()

        if options is not None:
            self._config.read_dict(options)
            self._has_config = True
            with open(self._config_file, 'w') as f:
                self._config.write(f)
        else:
            self._has_config = len(self._config.read(self._config_file))

    @property
    def raw_store(self):
        if not self._store:
            self._store = Store(self._config)
        return self._store

    @property
    def tracker(self):
        if self._tracker is None:
            self._tracker = DatasetTracker(self._config)
        return self._tracker

    @property
    def database(self):
        if self._database is None:
            self._database = DatabaseHelper(self._config)
        return self._database

    def cleanup(self):
        if self._database is not None:
            self._database.close()


@shell(context_settings={'help_option_names': ['-h', '--help']}, prompt="fzo> ")
@click.version_option(version=__version__, message=prog_notice)
@click.option(
    '--config',
    '-c',
    type=click.Path(exists=False),
    default='{}/config'.format(click.get_app_dir('fuzzyonions')),
    help="Path to an alternative configuration file.",
)
@click.pass_context
def main(ctx, config):
    """Helper scripts for the FlyBase scRNAseq project."""

    logging.basicConfig(format="fzo: %(module)s: %(message)s", level=logging.INFO)

    if not '/' in config and not exists(config):
        config = '{}/{}'.format(click.get_app_dir('fuzzyonions'), config)

    context = FzoContext(config)
    ctx.obj = context

    if not context.has_config:
        ctx.invoke(conf)

    ctx.call_on_close(context.cleanup)


@main.command()
@click.pass_obj
def ipython(ctx):
    """Start an interactive Python shell."""

    store = ctx.raw_store
    embed()


@main.command()
@click.pass_obj
def conf(ctx):
    """Edit the configuration."""

    if ctx.has_config:
        click.termui.edit(filename=ctx.config_file)
        ctx.reset()
    else:
        home = getenv('HOME')
        store_dir = f'{home}/scRNAseq/raw'
        proformae_dir = f'{home}/SVN_folders/proformae'

        store_dir = click.termui.prompt("Raw data directory", default=store_dir)
        proformae_dir = click.termui.prompt(
            "Proformae directory", default=proformae_dir
        )

        defaults = {
            'store': {'production': store_dir},
            'curation': {
                'proformae': proformae_dir,
                'trackfile': '{}/track.json'.format(click.get_app_dir('fuzzyonions')),
            },
        }
        ctx.reset(options=defaults)


main.add_command(store)
main.add_command(curate)
main.add_command(explorer)
main.add_command(tracker)
main.add_command(discover)
main.add_command(annots)

if __name__ == '__main__':
    main()
