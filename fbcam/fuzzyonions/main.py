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

from enum import IntFlag
import logging
from os import getenv
from os.path import exists
from configparser import ConfigParser

import click
from click_shell import shell
from IPython import embed

from fbcam.fuzzyonions import __version__
from fbcam.fuzzyonions.scea import FileStore, CombinedFileStore
from fbcam.fuzzyonions.explorer import explorer
from fbcam.fuzzyonions.proformae import ProformaGeneratorBuilder
from fbcam.fuzzyonions.curation import CuratedDatasetFactory
from fbcam.fuzzyonions.tracker import tracker, DatasetTracker

prog_name = "fuzzyonions"
prog_notice = f"""\
{prog_name} {__version__}
Copyright © 2021 Damien Goutte-Gattat

This program is released under the terms of the 1-clause BSD licence.
"""


class SourceStore(IntFlag):
    BOTH = 0,
    PRODUCTION = 1,
    STAGING = 2


class FzoContext(object):

    def __init__(self, config_file, source):
        self._config_file = config_file
        self._config = ConfigParser()
        self._source = source

        self.reset()

    @property
    def has_config(self):
        return self._has_config

    @property
    def config_file(self):
        return self._config_file

    def reset(self, options=None):
        self._store = None
        self._dataset = None
        self._subset = None
        self._subset_filters = []
        self._curation_factory = None
        self._tracker = None

        self._config.clear()

        if options is not None:
            self._config.read_dict(options)
            self._has_config = True
            with open(self._config_file, 'w') as f:
                self._config.write(f)
        else:
            self._has_config = len(self._config.read(self._config_file));

    @property
    def raw_store(self):
        if not self._store:
            prod_dir = self._config.get('store', 'production', fallback=None)
            staging_dir = self._config.get('store', 'staging', fallback=None)

            if prod_dir and staging_dir and self._source == SourceStore.BOTH:
                self._store = CombinedFileStore(prod_dir, staging_dir)
            elif prod_dir and (self._source == SourceStore.BOTH or
                               self._source == SourceStore.PRODUCTION):
                self._store = FileStore(prod_dir)
            elif staging_dir and (self._source == SourceStore.BOTH or
                                  self._source == SourceStore.STAGING):
                self._store = FileStore(staging_dir, staging=True)
            else:
                raise Exception("Invalid store configuration.")
        return self._store

    @property
    def tracker(self):
        if self._tracker is None:
            track = self._config.get('tracking', 'file', fallback=None)
            if not track:
                track = '{}/track.json'.format(click.get_app_dir('fuzzyonions'))
            self._tracker = DatasetTracker(track)
        return self._tracker

    @property
    def dataset(self):
        if self._dataset is None:
            self._dataset = self.raw_store.datasets[0]
        return self._dataset

    @property
    def subset(self):
        if self._subset is None:
            return self.dataset.experiment_design
        else:
            return self._subset

    @subset.setter
    def subset(self, subset):
        self._subset = subset
        if subset is None:
            self._subset_filters.clear()

    def load_dataset(self, dsid):
        self._dataset = self.raw_store.get(dsid)

    @property
    def proformae_folder(self):
        return self._config.get('proformae', 'directory')

    @property
    def curation_factory(self):
        if self._curation_factory is None:
            self._curation_factory = CuratedDatasetFactory(self.raw_store)
        return self._curation_factory

    def filter_subset(self, column, value):
        self.subset = self.subset.loc[self.subset[column] == value]
        self._subset_filters.append([column, value])

    def get_filter_string(self):
        if len(self._subset_filters) == 0:
            return '(all)'
        else:
            return ' > '.join([b for _, b in self._subset_filters])


@shell(context_settings={'help_option_names': ['-h', '--help']},
       prompt="fzo> ")
@click.version_option(version=__version__, message=prog_notice)
@click.option('--config', '-c', type=click.Path(exists=False),
              default='{}/config'.format(click.get_app_dir('fuzzyonions')),
              help="Path to an alternative configuration file.")
@click.option('--production', '-p', is_flag=True, default=False,
              help="Only use data from production server.")
@click.option('--staging', '-s', is_flag=True, default=False,
              help="Only use data from staging server.")
@click.pass_context
def main(ctx, config, production, staging):
    """Helper scripts for the FlyBase scRNAseq project."""

    logging.basicConfig(format="fzo: %(module)s: %(message)s",
                        level=logging.INFO)

    if not '/' in config and not exists(config):
        config = '{}/{}'.format(click.get_app_dir('fuzzyonions'), config)

    if production and staging:
        raise click.ClickException("Cannot use both --production and "
                                   "--staging.")
    if production:
        source = SourceStore.PRODUCTION
    elif staging:
        source = SourceStore.STAGING
    else:
        source = SourceStore.BOTH

    context = FzoContext(config, source)
    ctx.obj = context

    if not context.has_config:
        ctx.invoke(conf)


@main.command()
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


@main.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--with-reads/--without-reads', default=True,
              help="Extract number of reads per biosample.")
@click.option('--text', '-t', is_flag=True, default=False,
              help="Produce a text output instead of JSON.")
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def extract(ctx, specfile, with_reads, text, output):
    """Extract curation data from a dataset.
    
    This command expects a JSON-formatted file describing how to extract
    from a dataset the bits of data that are required for FlyBase
    curation. It produces a similar JSON file enriched with the
    extracted values.
    
    \b
    Sample JSON input file:
    {
        "symbol": <symbol to use in FlyBase>
        "dataset_id": <SCEA dataset ID>,
        "cell_types_column": <name of the column for cell types>,
        "excluded_cell_types": <list of cell types to ignore>,
        "conditions": <list of columns used to assign cells to samples>,
        "samples": [
            {
                "symbol": <suffix to add to the dataset-level symbol>,
                "selectors": <list of values used to select cells>
            },
            <repeat for as many samples as needed>
        ]
    }
    """

    ds = ctx.curation_factory.from_specfile(specfile)
    ds.extract(with_reads)
    if text:
        ds.to_text(output)
    else:
        ds.to_json(output)


@main.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.option('--header', '-H', is_flag=True, default=False,
              help="Writes an uncommented header line.")
@click.pass_obj
def sumexpr(ctx, specfile, output, header):
    """Summarize expression data from a dataset.
    
    This command expects a JSON file similar to the one used by the
    'extract' command. It computes a per-cell type summarised gene
    expression matrix from a dataset.
    """

    ds = ctx.curation_factory.from_specfile(specfile)
    result = ds.summarise_expression()
    if not header:
        # Write a commented header line (needed for harvdev processing)
        output.write('#genes\t')
        output.write('\t'.join(result.columns))
        output.write('\n')
    result.to_csv(output, sep='\t', header=header)


@main.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def proforma(ctx, specfile, output):
    """Generate a proforma for a dataset."""

    builder = ProformaGeneratorBuilder(ctx.proformae_folder, output)
    ds = ctx.curation_factory.from_specfile(specfile)
    ds.to_proforma(builder)


@main.command()
@click.argument('specfile', type=click.File('r'))
@click.pass_obj
def fixscea(ctx, specfile):
    """Generate correction files for the SCEA."""

    ds = ctx.curation_factory.from_specfile(specfile)
    ds.generate_scea_files('experiment-design.with-fbids.tsv',
                           'celltypes-fbterms.tsv')


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

        store_dir = click.termui.prompt("Raw data directory",
                                        default=store_dir)
        proformae_dir = click.termui.prompt("Proformae directory",
                                            default=proformae_dir)

        defaults = {
            'store': {'production': store_dir},
            'proformae': {'directory': proformae_dir}
            }
        ctx.reset(options=defaults)


main.add_command(explorer)
main.add_command(tracker)

if __name__ == '__main__':
    main()
