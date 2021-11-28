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

import json
import logging
import numpy
import pandas


class CuratedDataset(object):
    """Represents a scRNAseq dataset with associated curation data."""

    def __init__(self, spec, dataset):
        """Creates a new instance.
        
        :param spec: a structure containing the curation data
        :param dataset: the actual dataset
        """

        self._spec = spec
        self._ds = dataset

    @property
    def cell_type_column(self):
        """Gets the name of the column with cell type information."""

        return self._spec.get('Cell types column', None)

    @property
    def excluded_cell_types(self):
        """Gets a list of cell types to exclude.
        
        Any cell type listed will be excluded from the generated
        proforma and summarised expression table.
        """

        return self._spec.get('Excluded cell types', [])

    @property
    def extracted(self):
        """Indicates whether data extraction has been performed."""

        return self._spec.get('Extracted', False)

    @extracted.setter
    def set_extracted(self, value):
        """Sets the extraction marker."""

        self._spec['Extracted'] = value

    def extract(self, with_reads=True):
        """Extracts data from raw data files.
        
        For each sample, this method will extract the number of cells in
        the sample, and the number of reads with the with_reads
        parameter is True. If cell type information is available, it
        will also extract the number of cells for each cell type.
        
        :param with_reads: if True (the default), extract the number
            of reads
        """

        if 'Corrections' in self._spec:
            self._ds.apply_corrections(self._spec['Corrections'])

        for sample in self._spec['Samples']:
            subset = self._get_sample_subset(sample)

            # The number of cells is simply the number of rows
            n_cells = len(subset)
            sample['Cells'] = n_cells

            # Same, but per cell type
            sample['Cell types'] = {}
            if self.cell_type_column is not None:
                for cell_type in subset[self.cell_type_column].unique():
                    if cell_type not in self.excluded_cell_types:
                        n = len(subset.loc[subset[self.cell_type_column] == cell_type])
                        if n > 0:
                            sample['Cell types'][cell_type] = n

            # Get the number of reads from the raw expression matrix
            if with_reads:
                mm = self._ds.raw_expression

                # In some datasets (at least E-GEOD-100058... for now), there
                # is a discrepancy between the experiment design table and the
                # expression matrix, where not all cell identifiers from the
                # experiment design table have a corresponding column in the
                # expression matrix. We need to remove those offending cell
                # IDs before we look up the number of reads. Such a
                # discrepancy MAY be an indicator that there's something fishy
                # with the dataset, so if we detect it, we print a warning
                # giving the extent of the discrepancy (how many cells are
                # missing).
                present = subset.loc[subset['Assay'].isin(mm.columns)]['Assay']
                diff = n_cells - len(present)
                if diff > 0:
                    logging.warn(f"{sample['Symbol']}: {diff}/{n_cells} cells "
                                  "were removed because they are absent from "
                                  "the expression matrix.")

                n_reads = mm.loc[:, present].sum().sum()
                sample['Reads'] = int(n_reads)

        self.extracted = True

    def summarise_expression(self):
        """Generates a summarised expression table."""

        if 'Corrections' in self._spec:
            self._ds.apply_corrections(self._spec['Corrections'])

        # Get the normalized expression in exploitable form
        # HACK: The matrix read by SciPy's mmread function is filled with
        # zeros. To replace them with NaNs, we need to transform the spare
        # matrix into a dense matrix. This is probably not very efficient,
        # but it seems good enough even with some of the largest datasets
        # currently available on SCEA.
        matrix = self._ds.normalised_expression.transpose()
        matrix = matrix.sparse.to_dense()
        matrix.replace(0.0, numpy.nan, inplace=True)

        # Join the expression matrix with the experiment design table to
        # associate cell IDs with cell types
        expd = self._ds.experiment_design.set_index('Assay')
        matrix = matrix.join(expd[self.cell_type_column], on='cells')

        result = None
        for sample in self._spec['Samples']:
            subset = self._get_sample_subset(sample)

            # Subset of the expression matrix for this sample
            sm = matrix.loc[subset['Assay']]

            # Loop through cell types in this sample
            cell_types = subset.loc[:, self.cell_type_column].dropna().unique()
            for cell_type in [c for c in cell_types if c not in self.excluded_cell_types]:
                # Subset of the expression matrix for this cell type
                smc = sm.loc[sm[self.cell_type_column] == cell_type,:].set_index(self.cell_type_column)

                # Mean expression
                means = smc.mean()

                # "Spread" of expression
                spreads = smc.count() / len(smc)

                # Build the result dataframe
                d = pandas.DataFrame(data={'mean_expr': means, 'spread': spreads})
                d['celltype'] = cell_type
                d['sample'] = self._get_cluster_symbol(sample, cell_type)
                if result is not None:
                    result = result.append(d)
                else:
                    result = d

        result.index.rename('genes', inplace=True)
        return result.dropna()

    def to_text(self, output):
        """Writes a text description of the dataset.
        
        :param output: a file-like object to write into
        """

        for sample in self._spec['Samples']:
            symbol = self._get_sample_symbol(sample)
            output.write(f"Sample {symbol}\n")
            output.write(f"  Cells: {sample['Cells']}\n")
            if 'Reads' in sample:
                output.write(f"  Reads: {sample['Reads']}\n")
            for cell_type, cells in sample['Cell types'].items():
                output.write(f"    {cell_type}: {cells}\n")

    def to_json(self, output):
        """Writes a JSON description of the dataset.
        
        :param output: a file-like object to write into
        """

        json.dump(self._spec, output, indent=2)

    def to_proforma(self, builder):
        """Generates a proforma for the dataset.
        
        :param builder: a ProformaGeneratorBuilder object
        """

        if not self.extracted:
            self.extract()

        generator = builder.get_generator(template='pub_mini')
        generator.fill_template()

        generator = builder.get_generator(template='dataset/project')
        fills = {
            'LC1a': self._spec['Symbol'],
            'LC2b': 'transcriptome ; FBcv:0003034',
            'LC99a': self._spec['Dataset ID'],
            'LC99b': 'EMBL-EBI Single Cell Expression Atlas Datasets'
            }
        generator.fill_template(fills)

        for sample in self._spec['Samples']:
            symbol = self._get_sample_symbol(sample)
            stage = sample['Stage']
            title = sample['Title']

            generator = builder.get_generator(template='dataset/biosample')
            fills = {
                'LC1a': symbol,
                'LC6g': title,
                'LC2b': 'isolated cells ; FBcv:0003047',
                'LC3': self._spec['Symbol'],
                'LC4g': f'<e><t>{stage}<a><s><note>',
                'LC11m': 'multi-individual sample ; FBcv:0003141\n' +
                         'cell isolation ; FBcv:0003170'
                }
            generator.fill_template(fills)

            generator = builder.get_generator(template='dataset/assay')
            fills = {
                'LC1a': symbol + '_seq',
                'LC6g': f'Single-cell RNA-seq of {title}',
                'LC2b': 'single-cell RNA-Seq ; FBcv:0009000',
                'LC3': self._spec['Symbol'],
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
                'LC3': self._spec['Symbol'],
                'LC14b': symbol + '_seq',
                'LC6e': sample['Cells'],
                'LC6f': 'Number of cells in sample retained for analysis'
                }
            generator.fill_template(fills)

            for cell_type, n in sample['Cell types'].items():
                generator = builder.get_generator(template='dataset/subresult')
                fills = {
                    'LC1a': self._get_cluster_symbol(sample, cell_type),
                    'LC6g': f'Clustering analysis of {title}, {cell_type}s cluster',
                    'LC2b': 'transcriptional cell cluster ; FBcv:0009003',
                    'LC3': symbol + '_seq_clustering',
                    'LC4g': f'<e><t>{stage}<a>{cell_type}<s><note>',
                    'LC6e': n,
                    'LC6f': 'Number of cells in cluster'
                    }
                generator.fill_template(fills)

        generator.write_terminator()

    def generate_scea_files(self, new_expd_file, ct_file):
        """Generates the files with SCEA corrections.
        
        :param new_expd_file: the name of the file the enriched experiment
            design table should be written to
        :param ct_file: the name of the file the cell type corrections
            should be written to
        """

        if not 'Corrections' in self._spec:
            logging.warn("No corrections found in spec")
            return

        self._ds.apply_corrections(self._spec['Corrections'], only_new=True, target='scea')
        self._ds.experiment_design.to_csv(new_expd_file, sep='\t')

        cell_type_column = self.cell_type_column
        if cell_type_column is None:
            cell_type_column = self._spec.get('Source cell types column', None)
            if cell_type_column is None:
                logging.warn("No cell type information to correct")
                return

        for correction in self._spec['Corrections']:
            if correction['Source'] != cell_type_column:
                continue

            if 'Target' in correction and correction['Target'] != 'scea':
                continue

            with open(ct_file, 'w') as f:
                f.write('Original term\tProposed new term\tComment\n')
                for old, new, comment in correction['Values']:
                    f.write(f'{old}\t{new}\t{comment}\n')

    def _get_sample_symbol(self, sample):
        return self._spec['Symbol'] + sample['Symbol']

    def _get_cluster_symbol(self, sample, cell_type):
        replace_rules = [
            ('embryonic/larval ', ''),
            ('embryonic ', ''),
            ('larval ', ''),
            ('adult ', ''),
            (' ', '_')
            ]
        for rule in replace_rules:
            cell_type = cell_type.replace(rule[0], rule[1])
        sample_symbol = self._get_sample_symbol(sample)
        return f'{sample_symbol}_seq_clustering_{cell_type}s'

    def _get_sample_subset(self, sample):
        columns = self._spec.get('Conditions', None)
        subset = self._ds.experiment_design
        if columns:
            selectors = sample['Selectors']
            for i in range(len(columns)):
                subset = subset.loc[subset[columns[i]] == selectors[i]]

        return subset


class CuratedDatasetFactory(object):

    def __init__(self, store):
        self._store = store

    def from_specfile(self, specfile):
        spec = json.load(specfile)
        dataset = self._store.get(spec['Dataset ID'])
        return CuratedDataset(spec, dataset)
