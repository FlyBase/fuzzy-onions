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
import click
from click_shell.core import make_click_shell

from fbcam.fuzzyonions.proformae import ProformaGeneratorBuilder, ProformaType


class SourceDataset(object):
    """A SCEA dataset that is used as a source for a FlyBase dataset."""

    def __init__(self, spec):
        self._accession = spec['accession']
        self._ct_column = spec.get('cell_types_column', None)
        self._input_ct_column = spec.get('input_cell_types_column', None)
        self._simple_ct = spec.get('simplified_cell_types', None)
        self._corrections = []
        for corrset in spec.get('corrections', []):
            self._corrections.append(CorrectionSet(corrset))
        self._data = None

    @property
    def accession(self):
        """The SCEA accession number for this dataset."""

        return self._accession

    @property
    def cell_type_column(self):
        """The name of the column containing output cell types."""

        return self._ct_column

    @property
    def has_cell_types(self):
        """Indicate whether cell type annotations are available."""

        return self._ct_column is not None

    @property
    def input_cell_type_column(self):
        """The name of the column containing input cell types."""

        return self._input_ct_column

    @property
    def data(self):
        """The object representing the SCEA raw data files."""

        return self._data

    def get_simplified_cell_type(self, cell_type):
        """Get a simplified cell type name.
        
        Given an original cell type name, this method returns a
        simplified name that may be used when generating dataset
        symbols.
        """

        if self._simple_ct is not None:
            return self._simple_ct.get(cell_type, cell_type)
        return cell_type

    def get_corrections(self, target):
        """Get all correction sets applicable to the specified target."""

        return [c for c in self._corrections if c.valid_for(target)]

    def apply_corrections(self, only_new=False, target='internal'):
        """Apply correction sets for the specified target.
        
        :param only_new: if True, only apply corrections that do not
            modify existing columns
        :param: target: intended target for the corrections
        :return: the number of applied correction sets
        """

        expd = self.data.experiment_design
        n = 0
        for corrset in self.get_corrections(target):
            dest = corrset.destination
            if dest is not None:
                # We are adding a new column
                expd[corrset.destination] = pandas.Series(dtype='string')
            else:
                # We modify an existing column
                if only_new:
                    continue
                dest = corrset.source

            for old, new, _ in corrset.values:
                if len(new) == 0:
                    new = pandas.NA
                expd.loc[expd[corrset.source] == old, dest] = new

            n += 1

        return n

    def feed_data(self, store):
        """Inject raw data from a store.
        
        This method uses the provided store to fetch the raw data
        for the current dataset.
        """

        self._data = store.get(self.accession)


class CorrectionSet(object):
    """A set of corrections to apply to the ExperimentDesign table."""

    def __init__(self, data):
        self._target = data.get('target', 'both')
        self._source = data['source']
        self._destination = data.get('destination', None)
        self._values = data['values']

    @property
    def source(self):
        """The name of the column to correct."""

        return self._source

    @property
    def destination(self):
        """The name of the column to write corrected data to.
        
        If None, the corrections are intended to be applied directly
        to the source column.
        """

        return self._destination

    @property
    def values(self):
        """The corrections proper.
        
        This is a list of tuples, where each tuple contains:
        - the value to correct in the source column;
        - the corrected value to write to the destination column;
        - an optional comment.
        """

        return self._values

    def valid_for(self, target):
        """Indicate whether the set is applicable to a target."""

        return self._target in [target, 'both']


class FlybaseReference(object):

    def __init__(self, spec):
        self._fbrf = spec.get('fbrf', None)
        self._journal = spec.get('journal', None)

    @property
    def fbrf(self):
        return self._fbrf

    @property
    def journal(self):
        return self._journal


class SourceLab(object):

    def __init__(self, spec):
        self._name = spec.get('name', None)
        self._url = spec.get('url', None)

    @property
    def name(self):
        return self._name

    @property
    def url(self):
        return self._url


class DatasetBase(object):

    def __init__(self, spec):
        self._symbol = spec['symbol']
        self._title = spec.get('title', None)
        self._desc = spec.get('description', None)

        cv_terms = spec.get('cv_terms', {})
        self._go_cc = cv_terms.get('go_cc', [])
        self._go_mf = cv_terms.get('go_mf', [])
        self._go_bp = cv_terms.get('go_bp', [])
        self._so = cv_terms.get('so', [])
        self._fbcv = cv_terms.get('fbcv', [])

        protocols = spec.get('protocols', {})
        self._prot_collection = protocols.get('collection', None)
        self._prot_preparation = protocols.get('preparation', None)
        self._prot_assay = protocols.get('assay', None)
        self._prot_analysis = protocols.get('analysis', None)

        self._species = spec.get('species', {}).get('main', 'Dmel')

        self._excluded_ct = spec.get('excluded_cell_types', [])

    @property
    def symbol(self):
        return self._symbol

    @property
    def title(self):
        return self._title

    @property
    def description(self):
        if hasattr(self, '_desc'):
            return self._desc
        else:
            return None

    @property
    def entity_type(self):
        pass

    @property
    def data_type(self):
        pass

    @property
    def go_cellular_components(self):
        return self._go_cc

    @property
    def go_molecular_functions(self):
        return self._go_mf

    @property
    def go_biological_processes(self):
        return self._go_bp

    @property
    def fbcv(self):
        return self._fbcv

    @property
    def collection_protocol(self):
        return self._prot_collection

    @property
    def preparation_protocol(self):
        return self._prot_preparation

    @property
    def assay_protocol(self):
        return self._prot_assay

    @property
    def analysis_protocol(self):
        return self._prot_analysis

    @property
    def species(self):
        return self._species

    @property
    def other_species(self):
        return []

    @property
    def has_count(self):
        return False

    @property
    def count_label(self):
        return None

    def to_proforma(self, g):
        self.make_proforma_header(g)
        g.write_field('LC13a', '\n'.join(self.go_cellular_components))
        g.write_field('LC13b', '\n'.join(self.go_molecular_functions))
        g.write_field('LC13c', '\n'.join(self.go_biological_processes))
        g.write_field('LC13d', '\n'.join(self._so))
        g.write_field('LC4a', self.species)
        g.write_field('LC4i', '\n'.join(self.other_species))
        g.write_field('LC6d', 'N')
        if self.has_count:
            g.write_field('LC6e', self.count)
            g.write_field('LC6f', self.count_label)
        g.write_field('LC11m', '\n'.join(self.fbcv))
        g.write_field('LC11a', self.collection_protocol)
        g.write_field('LC6b', self.preparation_protocol)
        g.write_field('LC11c', self.assay_protocol)
        g.write_field('LC11e', self.analysis_protocol)

    def make_proforma_header(self, g):
        g.write_header()
        g.write_field('LC1f', 'new')
        g.write_field('LC1a', self.symbol)
        g.write_field('LC6g', self.title[0].upper() + self.title[1:])
        g.write_field('LC6a', self.description)
        g.write_field('LC2a', self.entity_type)
        g.write_field('LC2b', self.data_type)


class Project(DatasetBase):

    def __init__(self, spec, min_cluster_size=0):
        super().__init__(spec)
        self._min_cluster_size = min_cluster_size

        self._reference = FlybaseReference(spec.get('reference', {}))
        self._lab = SourceLab(spec.get('creator', {}))

        self._sources = []

        self._samples = []
        for sample in spec.get('samples', []):
            self._samples.append(Biosample(sample, self))

        self._subprojects = []
        for subproject in spec.get('subprojects', []):
            self._subprojects.append(Subproject(subproject, self))

    @property
    def entity_type(self):
        return 'project ; FBcv:0003023'

    @property
    def data_type(self):
        return 'transcriptome ; FBcv:0003034'

    @property
    def sources(self):
        return self._sources

    @property
    def other_species(self):
        species = []
        for sample in self.get_all_samples():
            species.extend(sample.other_species)
        return species

    @property
    def subprojects(self):
        return self._subprojects

    @property
    def has_subprojects(self):
        return len(self._subprojects) > 0

    @property
    def samples(self):
        return self._samples

    def has_samples(self):
        return len(self._samples) > 0

    @property
    def reference(self):
        return self._reference

    @property
    def lab(self):
        return self._lab

    @property
    def min_cluster_size(self):
        return self._min_cluster_size

    def get_all_samples(self):
        if self.has_samples:
            return self.samples
        else:
            samples = []
            for subproject in self.subprojects:
                samples.extend(subproject.get_all_samples())
            return samples

    def extract_reads(self):
        sample_by_cell_id = {}
        all_samples = self.get_all_samples()
        for sample in all_samples:
            for cell_id in sample.subset['Assay']:
                if not cell_id in sample_by_cell_id:
                    sample_by_cell_id[cell_id] = [sample]
                else:
                    sample_by_cell_id[cell_id].append(sample)

        for source in self.sources:
            with source.data.get_expression_matrix(raw=True) as mmf:
                mmf.set_progress_callback(lambda p: logging.info(f"Reading expression matrix: {p}% complete..."))
                for _, cell, value in mmf:
                    samples = sample_by_cell_id.get(cell)
                    if samples:
                        for sample in samples:
                            sample.assay.count += value

        for sample in all_samples:
            sample.assay.count = int(sample.assay.count)

    def summarise_expression(self):
        # Pre-processing
        # We prepare a data structure for each cluster in each sample.
        # That structure will be filled as we read the expression table
        # below. This allows to read the table only once.

        clusters_by_cell_id = {}
        samples = self.get_all_samples()
        for sample in samples:
            for cluster in sample.assay.result.clusters:
                for cell_id in cluster.subset['Assay']:
                    if not cell_id in clusters_by_cell_id:
                        clusters_by_cell_id[cell_id] = [cluster]
                    else:
                        clusters_by_cell_id[cell_id].append(cluster)

        logging.info("Preprocessing complete.")

        # Processing of the expression table

        for source in self.sources:
            with source.data.get_expression_matrix(raw=False) as mmf:
                mmf.set_progress_callback(lambda p: logging.info(f"Reading expression table: {p}% complete..."))
                for fbgn, cell, value in mmf:

                    if cell not in clusters_by_cell_id:
                        continue

                    for cluster in clusters_by_cell_id[cell]:
                        cluster.expression[fbgn] = cluster.expression.get(fbgn, 0) + value
                        cluster.presence[fbgn] = cluster.presence.get(fbgn, 0) + 1

        logging.info("Processing complete.")

        # Post-processing
        # We have all the values we need in the cluster structures.
        # Now we just need to assemble a DataFrame with them.

        result = None
        for sample in samples:
            for cluster in sample.assay.result.clusters:
                d = pandas.DataFrame(data={'expression': pandas.Series(data=cluster.expression),
                                           'presence': pandas.Series(data=cluster.presence)
                                           })
                # We compute mean expression for each gene...
                d['mean_expr'] = d['expression'] / d['presence']
                # and 'extent of expression'
                d['spread'] = d['presence'] / cluster.count
                d['celltype'] = cluster.cell_type
                d['sample'] = cluster.symbol

                d = d.drop(columns=['expression', 'presence'])
                if result is not None:
                    result = result.append(d)
                else:
                    result = d

        logging.info("Post-processing complete.")

        return result

    def to_proforma(self, g):
        super().to_proforma(g)

        g.write_field('LC7a', 'The EMBL-EBI\'s Single Cell Expression Atlas provides cell-level annotations, clustering data, raw and normalised read counts, and putative marker genes.')
        for source in self.sources:
            g.write_field('LC99a', source.accession)
            g.write_field('LC99b', 'EMBL-EBI Single Cell Expression Atlas Datasets')
        g.write_field('LC8c', f'[{self.lab.name}]({self.lab.url})')
        g.write_separator()

        for subproject in self.subprojects:
            subproject.to_proforma(g)

        for sample in self.samples:
            sample.to_proforma(g)

        g.write_terminator()


class Subproject(DatasetBase):

    def __init__(self, spec, parent):
        super().__init__(spec)
        self._parent = parent

        self._samples = []
        for sample in spec.get('samples', []):
            self._samples.append(Biosample(sample, self))

        self._subprojects = []
        for subproject in spec.get('subprojects', []):
            self._subprojects.append(Subproject(subproject, self))

    @property
    def symbol(self):
        return self.project.symbol + '_' + self._symbol

    @property
    def entity_type(self):
        return 'project ; FBcv:0003023'

    @property
    def data_type(self):
        return 'transcriptome ; FBcv:0003034'

    @property
    def project(self):
        return self._parent

    @property
    def subprojects(self):
        return self._subprojects

    @property
    def has_subprojects(self):
        return len(self._subprojects) > 0

    @property
    def samples(self):
        return self._samples

    @property
    def has_samples(self):
        return len(self._samples) > 0

    @property
    def other_species(self):
        species = []
        for sample in self.get_all_samples():
            species.extend(sample.other_species)
        return species

    def get_all_samples(self):
        if self.has_samples:
            return self.samples
        else:
            samples = []
            for subproject in self.subprojects:
                samples.extend(subproject.get_all_samples())
            return samples

    def to_proforma(self, g):
        self.make_proforma_header(g)
        g.write_field('LC3b', self.project.symbol)

        for subproject in self.subprojects:
            subproject.to_proforma(g)

        for sample in self.samples:
            sample.to_proforma(g)


class Biosample(DatasetBase):

    def __init__(self, spec, parent):
        super().__init__(spec)
        self._project = parent

        self._is_sn = spec.get('is_single_nucleus', 'no') == 'yes'

        self._other_species = spec.get('species', {}).get('other', [])

        self._strain = spec.get('strain', None)
        self._genotype = spec.get('genotype', None)

        self._fbbt = spec.get('anatomical_part', None)
        self._fbdv = spec.get('stage', None)
        self._sex = spec.get('sex', None)

        self._fbcv = spec.get('cv_terms', []).get('fbcv_sample', [])

        self._conditions = spec.get('conditions')
        self._selectors = spec.get('selectors')
        self._source_id = spec.get('source', None)
        self._source = None

        self._assay = Assay(spec, self)
        self._subset = None

    @property
    def symbol(self):
        return self.project.symbol + '_' + self._symbol

    @property
    def entity_type(self):
        return 'biosample ; FBcv:0003024'

    @property
    def data_type(self):
        if self._is_sn:
            return 'isolated nuclei ; FBcv:0009004'
        else:
            return 'isolated cells ; FBcv:0003047'

    @property
    def project(self):
        return self._project

    @property
    def top_project(self):
        parent = self.project
        while hasattr(parent, 'project'):
            parent = self.project
        return parent

    @property
    def is_single_nucleus(self):
        return self._is_sn

    @property
    def other_species(self):
        return self._other_species

    @property
    def strain(self):
        return self._strain

    @property
    def genotype(self):
        return self._genotype

    @property
    def anatomical_part(self):
        return self._fbbt

    @property
    def developmental_stage(self):
        return self._fbdv

    @property
    def sex(self):
        return self._sex

    @property
    def assay(self):
        return self._assay

    @property
    def subset(self):
        if self._subset is None:
            self._subset = self._get_subset()
        return self._subset

    @property
    def source(self):
        if self._source is None:
            self._source = self._get_source()
        return self._source

    def exclude_cell_type(self, cell_type):
        p = self
        while True:
            if cell_type in p._excluded_ct:
                return True
            if hasattr(p, 'project'):
                p = p.project
            else:
                return False

    def to_proforma(self, g):
        self.make_proforma_header(g)

        g.write_field('LC3', self.project.symbol)
        g.write_field('LC4a', self.species)
        g.write_field('LC4i', '\n'.join(self.other_species))
        g.write_field('LC4h', self.strain)
        g.write_field('LC4f', self.genotype)
        g.write_field('LC4g', self._get_tap_statement())
        # TODO: LC12a, LC12b
        g.write_field('LC6d', 'N')
        g.write_field('LC11m', '\n'.join(self.fbcv))
        g.write_field('LC11a', self.collection_protocol)

        g.write_separator()
        self.assay.to_proforma(g)

    def _get_tap_statement(self):
        t = self.developmental_stage
        if self.sex is not None:
            t += ' | ' + self.sex
        return f'<e><t>{t}<a>{self.anatomical_part}<s><note>'

    def _get_subset(self):
        subset = self.source.data.experiment_design
        if self._conditions:
            for i in range(len(self._conditions)):
                if isinstance(self._selectors[i], list):
                    subset = subset.loc[subset[self._conditions[i]].isin(self._selectors[i])]
                elif self._selectors[i] == '*':
                    pass
                else:
                    subset = subset.loc[subset[self._conditions[i]] == self._selectors[i]]

        return subset

    def _get_source(self):
        if self._source_id is None:
            return self.top_project.sources[0]
        else:
            return [s for s in self.top_project.sources if s.accession == self._source_id][0]


class Assay(DatasetBase):

    def __init__(self, spec, sample):
        if sample.is_single_nucleus:
            self._title = "Single-nucleus RNA-seq of " + sample.title
        else:
            self._title = "Single-cell RNA-seq of " + sample.title

        self._tech_ref = spec.get('technical_reference')
        self._biol_ref = spec.get('biological_reference')
        self._fbcv = spec.get('cv_terms', {}).get('fbcv_assay', [])

        self._sample = sample
        self._result = Result(self)
        self.count = 0

    @property
    def symbol(self):
        return self.sample.symbol + '_seq'

    @property
    def entity_type(self):
        return 'assay ; FBcv:0003025'

    @property
    def data_type(self):
        if self.sample.is_single_nucleus:
            return 'single-nucleus RNA-Seq ; FBcv:0009001'
        else:
            return 'single-cell RNA-Seq ; FBcv:0009000'

    @property
    def has_count(self):
        return True

    @property
    def count_label(self):
        return 'reads'

    @property
    def technical_reference(self):
        return self._tech_ref

    @property
    def biological_reference(self):
        return self._biol_ref

    @property
    def sample(self):
        return self._sample

    @property
    def result(self):
        return self._result

    def to_proforma(self, g):
        self.make_proforma_header(g)

        g.write_field('LC3', self.sample.project.symbol)
        g.write_field('LC14a', self.sample.symbol)
        # TODO: LC14d, LC14e
        g.write_field('LC4a', self.sample.species)
        g.write_field('LC6d', 'N')
        g.write_field('LC6e', self.count)
        g.write_field('LC6f', self.count_label)
        g.write_field('LC11m', '\n'.join(self.fbcv))

        g.write_separator()
        self.result.to_proforma(g)


class Result(DatasetBase):

    def __init__(self, assay):
        self._assay = assay
        self._title = "Clustering analysis of " + assay.sample.title
        self._clusters = None

    @property
    def symbol(self):
        return self.assay.symbol + '_clustering'

    @property
    def entity_type(self):
        return 'result ; FBcv:0003026'

    @property
    def data_type(self):
        return 'cell clustering analysis ; FBcv:0009002'

    @property
    def assay(self):
        return self._assay

    @property
    def count(self):
        return len(self.assay.sample.subset)

    @property
    def has_count(self):
        return True

    @property
    def count_label(self):
        return 'cells'

    @property
    def clusters(self):
        if self._clusters is None:
            self._clusters = self._get_clusters()
        return self._clusters

    def to_proforma(self, g):
        self.make_proforma_header(g)

        g.write_field('LC3', self.assay.sample.project.symbol)
        g.write_field('LC14b', self.assay.symbol)
        g.write_field('LC14h', 'Dmel R6.32')  # TODO: un-hardcode
        g.write_field('LC4a', self.assay.sample.species)
        g.write_field('LC6d', 'N')
        g.write_field('LC6e', self.count)
        g.write_field('LC6f', self.count_label)

        g.write_separator()
        for cluster in self.clusters:
            cluster.to_proforma(g)

    def _get_clusters(self):
        sample = self.assay.sample

        if not sample.source.has_cell_types:
            return None

        clusters = []
        for cell_type in sample.subset[sample.source.cell_type_column].dropna().unique():
            n = len(sample.subset.loc[sample.subset[sample.source.cell_type_column] == cell_type])
            if not self._exclude_cluster(cell_type, n):
                clusters.append(Cluster(cell_type, n, self))

        return clusters

    def _exclude_cluster(self, cell_type, count):
        if count < self.assay.sample.top_project.min_cluster_size:
            return True

        return self.assay.sample.exclude_cell_type(cell_type)


class Cluster(DatasetBase):

    def __init__(self, cell_type, count, result):
        self._cell_type = cell_type
        self._count = count
        self._result = result
        self._simple_ct = None
        self.expression = {}
        self.presence = {}

    @property
    def title(self):
        return f'{self.result.title}, {self.simplified_cell_type}s cluster'

    @property
    def entity_type(self):
        return 'result ; FBcv:0003026'

    @property
    def data_type(self):
        return 'transcriptional cell cluster ; FBcv:0009003'

    @property
    def result(self):
        return self._result

    @property
    def count(self):
        return self._count

    @property
    def has_count(self):
        return True

    @property
    def count_label(self):
        return 'cells'

    @property
    def cell_type(self):
        return self._cell_type

    @property
    def symbol(self):
        sct = self.simplified_cell_type
        rules = [
            (' ', '_'),
            ('/', '_'),
            (' ', '_')
            ]
        for rule in rules:
            sct = sct.replace(rule[0], rule[1])
        return f'{self.result.symbol}_{sct}s'

    @property
    def simplified_cell_type(self):
        if self._simple_ct is None:
            sct = self.result.assay.sample.source.get_simplified_cell_type(self._cell_type)
            rules = [
                ('embryonic/larval ', ''),
                ('embryonic ', ''),
                ('larval ', ''),
                ('adult ', '')
                ]
            for rule in rules:
                sct = sct.replace(rule[0], rule[1])
            self._simple_ct = sct
        return self._simple_ct

    @property
    def subset(self):
        sample = self.result.assay.sample
        ct_column = sample.source.cell_type_column
        cl_subset = sample.subset.loc[sample.subset[ct_column] == self.cell_type]
        return cl_subset

    def to_proforma(self, g):
        self.make_proforma_header(g)

        g.write_field('LC3', self.result.symbol)
        g.write_field('LC4a', self.result.assay.sample.species)
        g.write_field('LC4g', self._get_tap_statement())
        g.write_field('LC6d', 'N')
        g.write_field('LC6e', self.count)
        g.write_field('LC6f', self.count_label)

        g.write_separator()

    def _get_tap_statement(self):
        sample = self.result.assay.sample
        t = sample.developmental_stage
        if sample.sex is not None:
            t += ' | ' + sample.sex
        return f'<e><t>{t}<a>{self._cell_type}<s><note>'


class CuratedDataset(object):
    """Represents a scRNAseq dataset with associated curation data."""

    def __init__(self, spec, dataset, no_exclude=False, min_cluster=0):
        """Creates a new instance.
        
        :param spec: a structure containing the curation data
        :param dataset: the actual dataset
        """

        self._spec = spec
        self._ds = dataset
        self._exclude_cell_types = not no_exclude
        self._min_cluster = min_cluster

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

    def is_cell_type_excluded(self, cell_type, sample=None, size=-1):
        """Indicates whether a cell type should be excluded.

        A cell type can be excluded if it is specified in a
        `exclude_cell_types` directory, either directly in the dataset
        description object (where it applies to all samples in the
        dataset), or in the sample description object. Additionally,
        it can be excluded if the corresponding cluster contains less
        than a pre-defined threshold number of cells.
        """

        if size != -1 and size < self._min_cluster:
            return True

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
                    n = len(subset.loc[subset[self.cell_type_column] == cell_type])
                    if not self.is_cell_type_excluded(cell_type, sample, n) and n > 0:
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

            for cell_type in cell_types:
                simplified_cell_type = self._get_simplified_cell_type(cell_type)

                # Get the cells for this cluster in this sample
                ct_subset = subset.loc[subset[self.cell_type_column] == cell_type]
                if self.is_cell_type_excluded(cell_type, sample, len(ct_subset)):
                    continue

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
            stage = sample.get('stage', self._spec.get('all_samples', {}).get('stage', 'TODO: stage'))
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
            sample_terms = '\n'.join(sample.get('cv_terms', self._spec.get('all_samples', {}).get('cv_terms', ['<DEFAULT>'])))
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
            all_assays = self._spec.get('all_assays', {})
            assay = sample.get('assay', {})
            assay_terms = '\n'.join(assay.get('cv_terms', all_assays.get('cv_terms', ['<DEFAULT>'])))
            technical_ref = assay.get('technical_reference', all_assays.get('technical_reference', '<DEFAULT>'))
            if technical_ref != '' and technical_ref != '<DEFAULT>':
                technical_ref = self._spec['symbol'] + technical_ref + '_seq'
            biological_ref = assay.get('biological_reference', all_assays.get('biological_reference', '<DEFAULT>'))
            if biological_ref != '' and biological_ref != '<DEFAULT>':
                biological_ref = self._spec['symbol'] + biological_ref + '_seq'
            fills = {
                'LC1a': symbol + '_seq',
                'LC6g': assay_title,
                'LC2b': assay_type,
                'LC3': self._spec['symbol'],
                'LC14a': symbol,
                'LC6e': sample['reads'],
                'LC6f': 'reads',
                'LC14d': technical_ref,
                'LC14e': biological_ref,
                'LC11m': assay_terms
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


class CurationContext(object):

    def __init__(self, config, store, no_exclude, min_cluster_size):
        self._store = store
        self._no_exclude = no_exclude
        self._min_cluster = min_cluster_size
        self._proformae_dir = config.get('curation', 'proformae')

    def from_specfile(self, specfile):
        spec = json.load(specfile)
        dataset = self._store.get(spec['dataset_id'])
        return CuratedDataset(spec, dataset, self._no_exclude,
                              self._min_cluster)

    def get_proforma_builder(self, output):
        return ProformaGeneratorBuilder(self._proformae_dir, output)

    def source_dataset_from_specfile(self, specfile):
        if isinstance(specfile, str):
            with open(specfile, 'r') as f:
                spec = json.load(f)
        else:
            spec = json.load(specfile)
        dataset = SourceDataset(spec)
        dataset.feed_data(self._store)
        return dataset

    def dataset_from_specfile(self, specfile):
        spec = json.load(specfile)
        self._expand_defaults(spec)
        dataset = Project(spec, self._min_cluster)

        for source in spec['sources']:
            s = self.source_dataset_from_specfile(source)
            s.apply_corrections()
            dataset.sources.append(s)

        return dataset

    def _expand_defaults(self, spec, defaults={}):
        defaults = spec.get('defaults', defaults)

        for subproject in spec.get('subprojects', []):
            self._expand_defaults(subproject, defaults)

        for sample in spec.get('samples', []):
            self._copy_dict(defaults, sample)

    def _copy_dict(self, source, dest):
        for key, value in source.items():
            if key not in dest:
                dest[key] = value
            elif isinstance(value, dict):
                self._copy_dict(value, dest[key])


@click.group(invoke_without_command=True)
@click.option('--no-excluded-cell-types', '-n', 'no_exclude',
              is_flag=True, default=False,
              help="Do not exclude any cell types.")
@click.option('--min-cluster-size', '-m', default=0,
              help="Exclude cluster with less cells than specified.")
@click.pass_context
def curate(ctx, no_exclude, min_cluster_size):
    """Access the curation commands."""

    ctx.obj = CurationContext(ctx.obj.config, ctx.obj.raw_store, no_exclude,
                              min_cluster_size)
    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-curate> ")
        shell.cmdloop()


@curate.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--with-reads/--without-reads', default=True,
              help="Extract number of reads per biosample.")
@click.option('--text', '-t', is_flag=True, default=False,
              help="Produce a text output instead of JSON.")
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def extract(ctx, specfile, with_reads, text, output):
    """Extract curation data from a dataset."""

    ds = ctx.dataset_from_specfile(specfile)
    if with_reads:
        ds.extract_reads()

    if text:
        for sample in ds.get_all_samples():
            output.write(f"Sample {sample.symbol}\n")
            output.write(f"  Cells: {sample.assay.result.count}\n")
            if with_reads:
                output.write(f"  Reads: {sample.assay.count}\n")
            for cluster in sample.assay.result.clusters:
                output.write(f"    {cluster.cell_type}: {cluster.count}\n")
    else:
        # TODO: JSON output to be re-implemented
        pass


@curate.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.option('--header', '-H', is_flag=True, default=False,
              help="Writes an uncommented header line.")
@click.pass_obj
def sumexpr(ctx, specfile, output, header):
    """Summarize expression data from a dataset."""

    ds = ctx.dataset_from_specfile(specfile)
    result = ds.summarise_expression()
    if not header:
        # Write a commented header line (needed for harvdev processing)
        output.write('#genes\t')
        output.write('\t'.join(result.columns))
        output.write('\n')
    result.to_csv(output, sep='\t', header=header, float_format='%.6f')


@curate.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def proforma(ctx, specfile, output):
    """Generate a proforma for a dataset."""

    ds = ctx.dataset_from_specfile(specfile)
    ds.extract_reads()
    builder = ctx.get_proforma_builder(output)
    g = builder.get_generator(proforma_type=ProformaType.PUB_MINI)
    g.write_header()
    g.write_field('P22', ds.reference.fbrf)
    g.write_field('P2', ds.reference.journal)
    g.write_field('P19')
    g.write_separator()

    ds.to_proforma(builder.get_generator(proforma_type=ProformaType.DATASET))


@curate.command()
@click.argument('specfile', type=click.File('r'))
@click.pass_obj
def fixscea(ctx, specfile):
    """Generate correction files for the SCEA."""

    ds = ctx.source_dataset_from_specfile(specfile)
    new_expd_file = 'experiment-design.with-fbids.tsv'
    ct_file = 'celltypes-fbterms.tsv'

    if ds.apply_corrections(only_new=True, target='scea') > 0:
        ds.data.experiment_design.to_csv(new_expd_file, sep='\t')

        ct_column = ds.cell_type_column
        if ct_column is None:
            ct_column = ds.input_cell_type_column
            if ct_column is None:
                logging.warn("No cell type informations to correct")
                return

        for corrset in ds.get_corrections('scea'):
            if corrset.source != ct_column:
                continue

            with open(ct_file, 'w') as f:
                f.write('Original term\tProposed new term\tComment\n')
                for old, new, comment in corrset.values:
                    f.write(f'{old}\t{new}\t{comment}\n')

