# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2021 Damien Goutte-Gattat
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
from configparser import ConfigParser

import click
from click_shell import shell
from IPython import embed

from fbcam.fuzzyonions import __version__
from fbcam.fuzzyonions.scea import FileStore

prog_name = "fuzzyonions"
prog_notice = f"""\
{prog_name} {__version__}
Copyright © 2021 Damien Goutte-Gattat

This program is released under the terms of the 1-clause BSD licence.
"""


class FzoContext(object):

    def __init__(self, config_file):
        self._config = ConfigParser()
        self._config.read(config_file)

        self._store = None

    @property
    def raw_store(self):
        if not self._store:
            d = self._config.get('store', 'directory')
            self._store = FileStore(d)
        return self._store


@shell(context_settings={'help_option_names': ['-h', '--help']},
       prompt="fzo >")
@click.version_option(version=__version__, message=prog_notice)
@click.option('--config', '-c', type=click.Path(exists=True),
              default='{}/config'.format(click.get_app_dir('fuzzyonions')),
              help="Path to an alternative configuration file.")
@click.pass_context
def main(ctx, config):
    """Helper scripts for the FlyBase scRNAseq project."""

    logging.basicConfig(format="fzo: %(module)s: %(message)s",
                        level=logging.INFO)

    context = FzoContext(config)
    ctx.obj = context


@main.command()
@click.argument('dsid')
@click.pass_obj
def download(ctx, dsid):
    """Download a SCEA dataset.
    
    This command fetches the data for the specified dataset from the
    SCEA server.
    """

    ctx.raw_store.get(dsid)


@main.command('list')
@click.pass_obj
def list_datasets(ctx):
    """List the datasets in the local store."""

    n = len(ctx.raw_store.datasets)
    if n == 0:
        print("No local datasets available.")
        return

    print(f"{n} dataset(s) available:")
    for d in ctx.raw_store.datasets:
        print(f"{d.id}")


@main.command()
@click.pass_obj
def ipython(ctx):
    """Start an interactive Python shell."""

    store = ctx.raw_store
    embed()


if __name__ == '__main__':
    main()
