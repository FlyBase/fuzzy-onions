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

from fbcam.fuzzyonions.matrixmarket import MatrixMarketFile


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

        return self._spec.get('cell_types_column', None)

    @property
    def excluded_cell_types(self):
        """Gets a list of cell types to exclude.
        
        Any cell type listed will be excluded from the generated
        proforma and summarised expression table.
        """

        return self._spec.get('excluded_cell_types', [])

    @property
    def simplified_cell_types(self):
        """Gets a table of simplified cell type names.
        
        Gets a dictionary associating a FBbt cell type name to a
        simplified cell type name to be used when generating
        cluster symbols.
        """

        return self._spec.get('simplified_cell_types', {})

    @property
    def extracted(self):
        """Indicates whether data extraction has been performed."""

        return self._spec.get('extracted', False)

    @extracted.setter
    def extracted(self, value):
        """Sets the extraction marker."""

        self._spec['extracted'] = value

    def extract(self, with_reads=True):
        """Extracts data from raw data files.
        
        For each sample, this method will extract the number of cells in
        the sample, and the number of reads with the with_reads
        parameter is True. If cell type information is available, it
        will also extract the number of cells for each cell type.
        
        :param with_reads: if True (the default), extract the number
            of reads
        """

        if 'corrections' in self._spec:
            self._ds.apply_corrections(self._spec['corrections'])

        for sample in self._spec['samples']:
            subset = self._get_sample_subset(sample)

            # The number of cells is simply the number of rows
            n_cells = len(subset)
            sample['cells'] = n_cells

            # Same, but per cell type
            sample['cell_types'] = {}
            if self.cell_type_column is not None:
                for cell_type in subset[self.cell_type_column].unique():
                    if cell_type not in self.excluded_cell_types:
                        n = len(subset.loc[subset[self.cell_type_column] == cell_type])
                        if n > 0:
                            sample['cell_types'][cell_type] = n

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
                    logging.warn(f"{sample['symbol']}: {diff}/{n_cells} cells "
                                  "were removed because they are absent from "
                                  "the expression matrix.")

                n_reads = mm.loc[:, present].sum().sum()
                sample['reads'] = int(n_reads)

        self.extracted = True

    def summarise_expression(self):
        """Generates a summarised expression table.
        
        This method produces the summarised expression table by
        "manually" parsing the normalised expression table provided
        by the SCEA and doing the necessary computations on the fly.
        
        This is required for datasets with an expression table that
        is too big to fit entirely in memory.
        """

        if 'corrections' in self._spec:
            self._ds.apply_corrections(self._spec['corrections'])

        # Pre-processing
        # We prepare a data structure for each cluster in each sample.
        # That structure will be filled as we read the expression table
        # below. This allows to read the table only once.

        clusters = []
        clusters_by_cell_id = {}
        for sample in self._spec['samples']:
            subset = self._get_sample_subset(sample)

            # Get all cell types present in this sample
            cell_types = subset.loc[:, self.cell_type_column].dropna().unique()

            for cell_type in [c for c in cell_types if c not in self.excluded_cell_types]:
                simplified_cell_type = self._get_simplified_cell_type(cell_type)

                # Get the cells for this cluster in this sample
                ct_subset = subset.loc[subset[self.cell_type_column] == cell_type]

                # Fill the cluster structure
                cluster = {
                    'total_cells': len(ct_subset),
                    'symbol': self._get_cluster_symbol(sample, simplified_cell_type),
                    'cell_type': cell_type,
                    'expression': {},
                    'presence': {}
                    }
                clusters.append(cluster)

                # Allow for fast look up of which cluster a cell belongs to
                for cell_id in ct_subset['Assay']:
                    clusters_by_cell_id[cell_id] = cluster

        logging.info("Preprocessing complete.")

        # Processing of the expression table

        matrix_file = self._ds.get_expression_matrix_fullname(raw=False)
        with MatrixMarketFile(matrix_file) as mmf:
            mmf.set_progress_callback(lambda p: logging.info(f"Reading expression table: {p}% complete..."))
            for fbgn, cell, value in mmf:

                if cell not in clusters_by_cell_id:
                    continue

                cluster = clusters_by_cell_id[cell]
                cluster['expression'][fbgn] = cluster['expression'].get(fbgn, 0) + value
                cluster['presence'][fbgn] = cluster['presence'].get(fbgn, 0) + 1

        logging.info("Processing complete.")

        # Post-processing
        # We have all the values we need in the cluster structures.
        # Now we just need to assemble a DataFrame with them.

        result = None
        for cluster in clusters:
            n = cluster['total_cells']
            d = pandas.DataFrame(data={'expression': pandas.Series(data=cluster['expression']),
                                       'presence': pandas.Series(data=cluster['presence'])
                                       })
            # We compute mean expression for each gene...
            d['mean_expr'] = d['expression'] / d['presence']
            # and 'extent of expression'
            d['spread'] = d['presence'] / n
            d['cell_type'] = cluster['cell_type']
            d['symbol'] = cluster['symbol']

            d = d.drop(columns=['expression', 'presence'])

            if result is not None:
                result = result.append(d)
            else:
                result = d

        logging.info("Post-processing complete.")

        return result

    def to_text(self, output):
        """Writes a text description of the dataset.
        
        :param output: a file-like object to write into
        """

        for sample in self._spec['samples']:
            symbol = self._get_sample_symbol(sample)
            output.write(f"Sample {symbol}\n")
            output.write(f"  Cells: {sample['cells']}\n")
            if 'reads' in sample:
                output.write(f"  Reads: {sample['reads']}\n")
            for cell_type, cells in sample['cell_types'].items():
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

        with_cell_types = False
        for sample in self._spec['samples']:
            if len(sample['cell_types']):
                with_cell_types = True

        generator = builder.get_generator(template='pub_mini')
        fills = {
            'P22': 'TODO: FBrfXXXXXXX',
            'P2': 'TODO: <Journal abbreviation>'
        }
        generator.fill_template(fills)

        generator = builder.get_generator(template='dataset/project')
        fills = {
            'LC1a': self._spec['symbol'],
            'LC6g': 'TODO: Single-cell RNA-seq study of <tissue and stage> [upon condition]',
            'LC2b': 'transcriptome ; FBcv:0003034',
            'LC99a': self._spec['dataset_id'],
            'LC99b': 'EMBL-EBI Single Cell Expression Atlas Datasets',
            'LC7a': 'The EMBL-EBI\'s Single Cell Expression Atlas provides cell-level annotations, clustering data, raw and normalised read counts, and putative marker genes.',
            'LC8c': 'TODO: [Name of the research group](URL of their website)',
            'LC13c': 'TODO: [GO terms for biological processes studied]',
            'LC11m':' TODO: <FBcv terms from \'study design\' (FBcv:0003130)>',
            'LC6b': 'TODO: <How the biological material was processed>',
            'LC11c': 'TODO: <How the mRNA libraries were sequenced>',
            'LC11e': 'TODO: <How the sequencing data were analysed>'
            }
        if with_cell_types:
            fills['LC6a'] = 'TODO: A characterization of the diverse populations of cells in...'
        else:
            fills['LC6a'] = 'TODO: A single-cell transcriptomic study of the diverse populations from...'
        generator.fill_template(fills)

        for sample in self._spec['samples']:
            symbol = self._get_sample_symbol(sample)
            stage = sample['stage']
            title = sample['title']

            generator = builder.get_generator(template='dataset/biosample')
            fills = {
                'LC1a': symbol,
                'LC6g': title[0].upper() + title[1:],
                'LC2b': 'isolated cells ; FBcv:0003047',
                'LC3': self._spec['symbol'],
                'LC4g': f'<e><t>{stage}<a>TODO: tissue<s><note>',
                'LC11m': 'TODO: <FBcv terms from \'biosample attribute\' (FBcv:0003137)>',
                'LC4h': 'TODO: [as needed]',
                'LC4f': 'TODO: [as needed]',
                'LC12a': 'TODO: [as needed]',
                'LC12b': 'TODO: [as needed]',
                'LC11a': 'TODO: <How the sample was obtained>'
                }
            generator.fill_template(fills)

            generator = builder.get_generator(template='dataset/assay')
            fills = {
                'LC1a': symbol + '_seq',
                'LC6g': f'Single-cell RNA-seq of {title}',
                'LC2b': 'single-cell RNA-Seq ; FBcv:0009000',
                'LC3': self._spec['symbol'],
                'LC14a': symbol,
                'LC6e': sample['reads'],
                'LC6f': 'reads',
                'LC14d': 'TODO: [LC1a symbol of technical reference assay, if relevant]',
                'LC14e': 'TODO: [LC1a symbol of biological reference assay, if relevant]',
                'LC11m': 'TODO: <FBcv terms from \'assay method\' (FBcv:0003208)>'
                }
            generator.fill_template(fills)

            generator = builder.get_generator(template='dataset/result')
            fills = {
                'LC1a': symbol + '_seq_clustering',
                'LC6g': f'Clustering analysis of {title}',
                'LC2b': 'cell clustering analysis ; FBcv:0009002',
                'LC3': self._spec['symbol'],
                'LC14b': symbol + '_seq',
                'LC6e': sample['cells'],
                'LC6f': 'cells'
                }
            generator.fill_template(fills)

            for cell_type, n in sample['cell_types'].items():
                generator = builder.get_generator(template='dataset/subresult')
                simplified_cell_type = self._get_simplified_cell_type(cell_type)
                fills = {
                    'LC1a': self._get_cluster_symbol(sample, simplified_cell_type),
                    'LC6g': f'Clustering analysis of {title}, {simplified_cell_type}s cluster',
                    'LC2b': 'transcriptional cell cluster ; FBcv:0009003',
                    'LC3': symbol + '_seq_clustering',
                    'LC4g': f'<e><t>{stage}<a>{cell_type}<s><note>',
                    'LC6e': n,
                    'LC6f': 'cells'
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

        if not 'corrections' in self._spec:
            logging.warn("No corrections found in spec")
            return

        if self._ds.apply_corrections(self._spec['corrections'], only_new=True, target='scea') > 0:
            self._ds.experiment_design.to_csv(new_expd_file, sep='\t')

        cell_type_column = self.cell_type_column
        if cell_type_column is None:
            cell_type_column = self._spec.get('source_cell_types_column', None)
            if cell_type_column is None:
                logging.warn("No cell type information to correct")
                return

        for correction in self._spec['corrections']:
            if correction['source'] != cell_type_column:
                continue

            if 'target' in correction and correction['target'] != 'scea':
                continue

            with open(ct_file, 'w') as f:
                f.write('Original term\tProposed new term\tComment\n')
                for old, new, comment in correction['values']:
                    f.write(f'{old}\t{new}\t{comment}\n')

    def _get_sample_symbol(self, sample):
        return self._spec['symbol'] + sample['symbol']

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
        columns = self._spec.get('conditions', None)
        subset = self._ds.experiment_design
        if columns:
            selectors = sample['selectors']
            for i in range(len(columns)):
                subset = subset.loc[subset[columns[i]] == selectors[i]]

        return subset

    def _get_simplified_cell_type(self, cell_type):
        return self.simplified_cell_types.get(cell_type, cell_type)


class CuratedDatasetFactory(object):

    def __init__(self, store):
        self._store = store

    def from_specfile(self, specfile):
        spec = json.load(specfile)
        dataset = self._store.get(spec['dataset_id'])
        return CuratedDataset(spec, dataset)
