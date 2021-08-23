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
import json
from configparser import ConfigParser

import click
from click_shell import shell
from IPython import embed
import numpy
import pandas

from fbcam.fuzzyonions import __version__
from fbcam.fuzzyonions.scea import FileStore
from fbcam.fuzzyonions.explorer import explorer
from fbcam.fuzzyonions.proformae import ProformaGeneratorBuilder

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
        self._dataset = None
        self._subset = None

    @property
    def raw_store(self):
        if not self._store:
            d = self._config.get('store', 'directory')
            self._store = FileStore(d)
        return self._store

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

    def load_dataset(self, dsid):
        self._dataset = self.raw_store.get(dsid)

    @property
    def proformae_folder(self):
        return self._config.get('proformae', 'directory')


@shell(context_settings={'help_option_names': ['-h', '--help']},
       prompt="fzo> ")
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
        "Symbol": <symbol to use in FlyBase>
        "Dataset ID": <SCEA dataset ID>,
        "Cell types column": <name of the column for cell types>,
        "Excluded cell types": <list of cell types to ignore>,
        "Conditions": <list of columns used to assign cells to samples>,
        "Samples": [
            {
                "Symbol": <suffix to add to the dataset-level symbol>,
                "Selectors": <list of values used to select cells>
            },
            <repeat for as many samples as needed>
        ]
    }
    """

    spec = json.load(specfile)
    ds = ctx.raw_store.get(spec['Dataset ID'])
    cell_type_column = spec.get('Cell types column', None)
    excluded_cell_types = spec.get('Excluded cell types', [])
    columns = spec.get('Conditions', None)

    if 'Corrections' in spec:
        ds.apply_corrections(spec['Corrections'])

    for sample in spec['Samples']:

        # Get the subset of cells for this sample
        subset = ds.experiment_design
        if columns:
            selectors = sample['Selectors']
            for i in range(len(columns)):
                subset = subset.loc[subset[columns[i]] == selectors[i]]

        # The number of cells is simply the number of rows
        sample['Cells'] = len(subset)

        # Same, but per cell type
        sample['Cell types'] = {}
        if cell_type_column is not None:
            for cell_type in subset[cell_type_column].unique():
                if cell_type not in excluded_cell_types:
                    n = len(subset.loc[subset[cell_type_column] == cell_type])
                    if n > 0:
                        sample['Cell types'][cell_type] = n

        # Get the number of reads from the raw expression matrix
        if with_reads:
            mm = ds.raw_expression
            nreads = mm.loc[:, subset['Assay']].sum().sum()
            sample['Reads'] = int(nreads)

    if text:
        for sample in spec['Samples']:
            symbol = spec['Symbol'] + sample['Symbol']
            output.write(f"Sample {symbol}\n")
            output.write(f"  Cells: {sample['Cells']}\n")
            if with_reads:
                output.write(f"  Reads: {sample['Reads']}\n")
            for cell_type, cells in sample['Cell types'].items():
                output.write(f"    {cell_type}: {cells}\n")
    else:
        json.dump(spec, output, indent=2)


@main.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def sumexpr(ctx, specfile, output):
    """Summarize expression data from a dataset.
    
    This command expects a JSON file similar to the one used by the
    'extract' command. It computes a per-cell type summarised gene
    expression matrix from a dataset.
    """

    spec = json.load(specfile)
    ds = ctx.raw_store.get(spec['Dataset ID'])
    cell_type_column = spec['Cell types column']
    excluded_cell_types = spec.get('Excluded cell types', [])
    columns = spec.get('Conditions', None)

    if 'Corrections' in spec:
        ds.apply_corrections(spec['Corrections'])

    # Get the normalized expression data in exploitable form
    # HACK: The matrix read by SciPy's mmread function is filled with
    # zeros. To replace them with NaN, we need to transform the spare
    # matrix into a dense matrix. This is probably not very efficient,
    # but it seems good enough even with some of the largest datasets
    # currently available on SCEA.
    matrix = ds.normalised_expression.transpose()
    matrix = matrix.sparse.to_dense()
    matrix.replace(0.0, numpy.nan, inplace=True)

    # Join the expression matrix with the experiment design table to
    # associate cell IDs with cell types
    expd = ds.experiment_design.set_index('Assay')
    matrix = matrix.join(expd[cell_type_column], on='cells')

    result = None
    for sample in spec['Samples']:

        # Get the subset of cells for this sample
        subset = ds.experiment_design
        if columns:
            selectors = sample['Selectors']
            for i in range(len(columns)):
                subset = subset.loc[subset[columns[i]] == selectors[i]]

        # Subset of the expression matrix for this sample
        sm = matrix.loc[subset['Assay']]

        # Loop through cell types in this sample
        cell_types = subset.loc[:, cell_type_column].dropna().unique()
        for cell_type in [c for c in cell_types if c not in excluded_cell_types]:
            # Subset of the expression matrix for this cell type
            smc = sm.loc[sm[cell_type_column] == cell_type,:].set_index(cell_type_column)

            # Mean expression
            means = smc.mean()

            # "Spread" of expression
            spreads = smc.count() / len(smc)

            # Build the result dataframe
            d = pandas.DataFrame(data={'mean_expr': means, 'spread': spreads})
            d['celltype'] = cell_type
            d['sample'] = spec['Symbol'] + sample['Symbol']
            if result is not None:
                result = result.append(d)
            else:
                result = d

    result.index.rename('genes', inplace=True)
    result.dropna().to_csv(output, sep='\t')


@main.command()
@click.argument('spec', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def proforma(ctx, spec, output):
    """Generate a proforma for a dataset."""

    spec = json.load(spec)
    builder = ProformaGeneratorBuilder(ctx.proformae_folder, output)

    generator = builder.get_generator(template='pub_mini')
    generator.fill_template()

    generator = builder.get_generator(template='dataset/project')
    fills = {
        'LC1a': spec['Symbol'],
        'LC2b': 'transcriptome ; FBcv:0003034',
        'LC99a': spec['Dataset ID'],
        'LC99b': 'EMBL-EBI Single Cell Expression Atlas Datasets'
        }
    generator.fill_template(fills)

    for sample in spec['Samples']:
        symbol = spec['Symbol'] + sample['Symbol']
        stage = sample['Stage']
        title = sample['Title']

        generator = builder.get_generator(template='dataset/biosample')
        fills = {
            'LC1a': symbol,
            'LC6g': title,
            'LC2b': 'isolated cells ; FBcv:0003047',
            'LC3': spec['Symbol'],
            'LC4g': f'<e><t>{stage}<a><s><note>',
            'LC6e': sample['Cells'],
            'LC6f': 'Number of cells in sample',
            'LC11m': 'multi-individual sample ; FBcv:0003141\n' +
                     'cell isolation ; FBcv:0003170'
            }
        generator.fill_template(fills)

        generator = builder.get_generator(template='dataset/assay')
        fills = {
            'LC1a': symbol + '_seq',
            'LC6g': f'Single-cell RNA-seq of {title}',
            'LC2b': 'single-cell RNA-Seq ; FBcv:0009000',
            'LC3': spec['Symbol'],
            'LC14a': symbol,
            'LC6e': sample['Reads'],
            'LC6f': 'Number of reads'
            }
        generator.fill_template(fills)

        generator = builder.get_generator(template='dataset/result')
        fills = {
            'LC1a': symbol + '_seq_clustering',
            'LC6g': f'Clustering analysis of {title}',
            'LC2b': 'cell clustering analysis ; FBcv:0009002',
            'LC3': spec['Symbol'],
            'LC14b': symbol + '_seq',
            'LC6d': 'Y'
            }
        generator.fill_template(fills)

        for cell_type, n in sample['Cell types'].items():
            ct_symbol = cell_type.replace(' ', '_')
            generator = builder.get_generator(template='dataset/subresult')
            fills = {
                'LC1a': f'{symbol}_cluster_{ct_symbol}s',
                'LC6g': f'Clustering analysis of {title}, {cell_type}s cluster',
                'LC2b': 'transcriptional cell cluster ; FBcv:0009003',
                'LC3': symbol + '_seq_clustering',
                'LC4g': f'<e><t>{stage}<a>{cell_type}<s><note>',
                'LC6e': n,
                'LC6f': 'Number of cells in cluster'
                }
            generator.fill_template(fills)

    generator.write_terminator()


@main.command()
@click.argument('spec', type=click.File('r'))
@click.pass_obj
def fixscea(ctx, spec):
    """Generate correction files for the SCEA."""

    spec = json.load(spec)
    if not 'Corrections' in spec:
        print("No corrections found in spec file")
        return

    ds = ctx.raw_store.get(spec['Dataset ID'])
    ds.apply_corrections(spec['Corrections'], only_new=True)
    ds.experiment_design.to_csv('experiment-design.with-fbids.tsv', sep='\t')

    cell_type_column = spec['Cell types column']
    for correction in spec['Corrections']:
        if correction['Source'] != cell_type_column:
            continue

        with open('celltypes-fbterms.tsv', 'w') as f:
            f.write('Original term\tProposed new term\tComment\n')
            for old, new, comment in correction['Values']:
                f.write(f'{old}\t{new}\t{comment}\n')


main.add_command(explorer)

if __name__ == '__main__':
    main()
