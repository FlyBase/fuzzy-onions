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
import pandas


class CuratedDataset(object):
    """Represents a scRNAseq dataset with associated curation data."""

    def __init__(self, spec, dataset, no_exclude=False):
        """Creates a new instance.
        
        :param spec: a structure containing the curation data
        :param dataset: the actual dataset
        """

        self._spec = spec
        self._ds = dataset
        self._exclude_cell_types = not no_exclude

    @property
    def cell_type_column(self):
        """Gets the name of the column with cell type information."""

        return self._spec.get('cell_types_column', None)

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

    def is_cell_type_excluded(self, cell_type, sample=None):
        """Indicates whether a cell type should be excluded.

        A cell type can be excluded if it is specified in a
        `exclude_cell_types` directory, either directly in the dataset
        description object (where it applies to all samples in the
        dataset), or in the sample description object.
        """

        if not self._exclude_cell_types:
            return False

        if sample and 'excluded_cell_types' in sample:
            if cell_type in sample['excluded_cell_types']:
                return True

        if 'excluded_cell_types' in self._spec:
            if cell_type in self._spec['excluded_cell_types']:
                return True

        return False

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

        if with_reads:
            sample_by_cell_id = {}

        for sample in self._spec['samples']:
            subset = self._get_sample_subset(sample)

            # The number of cells is simply the number of rows
            n_cells = len(subset)
            sample['cells'] = n_cells

            # Same, but per cell type
            sample['cell_types'] = {}
            if self.cell_type_column is not None:
                for cell_type in subset[self.cell_type_column].unique():
                    if not self.is_cell_type_excluded(cell_type, sample):
                        n = len(subset.loc[subset[self.cell_type_column] == cell_type])
                        if n > 0:
                            sample['cell_types'][cell_type] = n

            if with_reads:
                # Allow for fast look up of which sample a cell belongs to
                for cell_id in subset['Assay']:
                    if not cell_id in sample_by_cell_id:
                        sample_by_cell_id[cell_id] = [sample]
                    else:
                        sample_by_cell_id[cell_id].append(sample)
                sample['reads'] = 0

        # Get the number of reads from the raw expression matrix
        if with_reads:
            with self._ds.get_expression_matrix(raw=True) as mmf:
                mmf.set_progress_callback(lambda p: logging.info(f"Reading expression matrix: {p}% complete..."))
                for _, cell, value in mmf:
                    samples = sample_by_cell_id.get(cell)
                    if samples:
                        for sample in samples:
                            sample['reads'] += value

            for sample in self._spec['samples']:
                sample['reads'] = int(sample['reads'])

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

            for cell_type in [c for c in cell_types if not self.is_cell_type_excluded(c, sample)]:
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
                    if not cell_id in clusters_by_cell_id:
                        clusters_by_cell_id[cell_id] = [cluster]
                    else:
                        clusters_by_cell_id[cell_id].append(cluster)

        logging.info("Preprocessing complete.")

        # Processing of the expression table

        with self._ds.get_expression_matrix(raw=False) as mmf:
            mmf.set_progress_callback(lambda p: logging.info(f"Reading expression table: {p}% complete..."))
            for fbgn, cell, value in mmf:

                if cell not in clusters_by_cell_id:
                    continue

                for cluster in clusters_by_cell_id[cell]:
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
            d['celltype'] = cluster['cell_type']
            d['sample'] = cluster['symbol']

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
        fills = {}
        if 'reference' in self._spec:
            fills['P22'] = self._spec['reference'].get('fbrf', '<DEFAULT>')
            fills['P2'] = self._spec['reference'].get('journal', '<DEFAULT>')
        generator.fill_template(fills)

        generator = builder.get_generator(template='dataset/project')
        fills = {
            'LC1a': self._spec['symbol'],
            'LC99a': self._spec['dataset_id']
            }
        if with_cell_types:
            fills['LC6a'] = 'TODO: A characterization of the diverse populations of cells in...'
        else:
            fills['LC6a'] = 'TODO: A single-cell transcriptomic study of the diverse populations from...'
        if 'project' in self._spec:
            project = self._spec['project']
            fills['LC6g'] = project.get('title', '<DEFAULT>')
            fills['LC6a'] = project.get('description')
            fills['LC13c'] = '\n'.join(project.get('go_terms', {}).get('biological_process', []))
            fills['LC11m'] = '\n'.join(project.get('protocols', {}).get('cv_terms', []))
            fills['LC6b'] = project.get('protocols', {}).get('preparation')
            fills['LC11c'] = project.get('protocols', {}).get('assay')
            fills['LC11e'] = project.get('protocols', {}).get('analysis')
            source = project['source']
            fills['LC8c'] = f"[{source['name']}]({source['url']})"
        generator.fill_template(fills)

        for sample in self._spec['samples']:
            symbol = self._get_sample_symbol(sample)
            stage = sample.get('stage', self._spec['all_samples']['stage'])
            if 'sex' in sample:
                stage += ' | ' + sample['sex']
            title = sample['title']
            tissue = sample.get('anatomical_part', 'TODO: tissue')
            genotype, strain = self._get_strain_and_genotype(sample)
            single_nuclei = sample.get('is_single_nuclei', self._spec.get('all_samples', {}).get('is_single_nuclei', 'no')) == 'yes'
            if single_nuclei:
                sample_type = 'isolated nuclei ; FBcv:0009004'
                assay_type = 'single-nucleus RNA-Seq ; FBcv:0009001'
                assay_title = f'Single-nucleus RNA-seq of {title}'
            else:
                sample_type = 'isolated cells ; FBcv:0003047'
                assay_type = 'single-cell RNA-Seq ; FBcv:0009000'
                assay_title = f'Single-cell RNA-Seq of {title}'

            generator = builder.get_generator(template='dataset/biosample')
            sample_terms = '\n'.join(sample.get('cv_terms', self._spec.get('all_samples', {}).get('cv_terms', '<DEFAULT>')))
            sample_prep = sample.get('isolation', self._spec.get('all_samples', {}).get('isolation', '<DEFAULT>'))
            entities = sample.get('entities', self._spec.get('all_samples', {}).get('entities'))
            if entities is None:
                entity_names = '<DEFAULT>'
                entity_types = '<DEFAULT>'
            else:
                entity_names = '\n! LC12a. Experimental entity :'.join([a[0] for a in entities])
                entity_types = '\n! LC12b. Type of experimental entity :'.join([a[1] for a in entities])
            fills = {
                'LC1a': symbol,
                'LC6g': title[0].upper() + title[1:],
                'LC2b': sample_type,
                'LC3': self._spec['symbol'],
                'LC4g': f'<e><t>{stage}<a>{tissue}<s><note>',
                'LC11m': sample_terms,
                'LC4h': strain,
                'LC4f': genotype,
                'LC12a': entity_names,
                'LC12b': entity_types,
                'LC11a': sample_prep
                }
            generator.fill_template(fills)

            generator = builder.get_generator(template='dataset/assay')
            fills = {
                'LC1a': symbol + '_seq',
                'LC6g': assay_title,
                'LC2b': assay_type,
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
            (' ', '_'),
            ('-', '_'),
            ('/', '_')
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
                if isinstance(selectors[i], list):
                    subset = subset.loc[subset[columns[i]].isin(selectors[i])]
                elif selectors[i] == '*':
                    pass
                else:
                    subset = subset.loc[subset[columns[i]] == selectors[i]]

        return subset

    def _get_simplified_cell_type(self, cell_type):
        return self.simplified_cell_types.get(cell_type, cell_type)

    def _get_strain_and_genotype(self, sample):
        genotype = sample.get('genotype')
        if genotype:
            return (genotype, '')

        strain = sample.get('strain')
        if strain:
            return ('', strain)

        if 'all_samples' in self._spec:
            genotype = self._spec['all_samples'].get('genotype')
            if genotype:
                return (genotype, '')
            strain = self._spec['all_samples'].get('strain')
            if strain:
                return ('', strain)

        return ('TODO: [as needed]', 'TODO: [as needed]')


class CuratedDatasetFactory(object):

    def __init__(self, store, no_exclude=False):
        self._store = store
        self._no_exclude = no_exclude

    def from_specfile(self, specfile):
        spec = json.load(specfile)
        dataset = self._store.get(spec['dataset_id'])
        return CuratedDataset(spec, dataset, self._no_exclude)
