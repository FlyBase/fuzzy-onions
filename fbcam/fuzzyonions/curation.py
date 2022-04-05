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
        self._go_cc = cv_terms.get('go_cc')
        self._go_mf = cv_terms.get('go_mf')
        self._go_bp = cv_terms.get('go_bp')
        self._so = cv_terms.get('so')
        self._fbcv = cv_terms.get('fbcv')

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
    def so_terms(self):
        return self._so

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


class ProjectContainer(DatasetBase):

    def __init__(self, spec, parent=None):
        super().__init__(spec)

        self._parent = parent

        self._subprojects = []
        for subproject in spec.get('subprojects', []):
            self._subprojects.append(ProjectContainer(subproject, self))

        self._samples = []
        for sample in spec.get('samples', []):
            self._samples.append(Biosample(sample, self))

    @property
    def subprojects(self):
        return self._subprojects

    @property
    def has_subprojects(self):
        return len(self.subprojects) > 0

    @property
    def samples(self):
        return self._samples

    @property
    def has_samples(self):
        return len(self.samples) > 0

    @property
    def project(self):
        return self._parent

    @property
    def is_top_project(self):
        return self._parent is None

    @property
    def top_project(self):
        p = self
        while not p.is_top_project:
            p = self.project
        return p

    @property
    def symbol(self):
        if self.is_top_project:
            return self._symbol
        else:
            return self.project.symbol + '_' + self._symbol

    @property
    def entity_type(self):
        return 'project ; FBcv:0003023'

    @property
    def data_type(self):
        return 'transcriptome ; FBcv:0003034'

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


class Project(ProjectContainer):

    def __init__(self, spec, min_cluster_size=0):
        super().__init__(spec)
        self._min_cluster_size = min_cluster_size

        self._reference = FlybaseReference(spec.get('reference', {}))
        self._lab = SourceLab(spec.get('creator', {}))

        self._sources = []

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
    def reference(self):
        return self._reference

    @property
    def lab(self):
        return self._lab

    @property
    def min_cluster_size(self):
        return self._min_cluster_size

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


class Biosample(DatasetBase):

    def __init__(self, spec, parent):
        super().__init__(spec)
        self._project = parent
        if self._desc is None:
            self._desc = ""

        self._is_sn = spec.get('is_single_nucleus', 'no') == 'yes'

        self._other_species = spec.get('species', {}).get('other', [])

        self._strain = spec.get('strain', None)
        self._genotype = spec.get('genotype', None)

        self._fbbt = spec.get('anatomical_part', None)
        self._fbdv = spec.get('stage', None)
        self._sex = spec.get('sex', None)

        self._fbcv = spec.get('cv_terms', {}).get('fbcv_sample')

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
        return self.project.top_project

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
            if p.project is not None:
                p = p.project
            else:
                return False

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

        self._desc = ""

        self._tech_ref = spec.get('technical_reference')
        self._biol_ref = spec.get('biological_reference')
        self._fbcv = spec.get('cv_terms', {}).get('fbcv_assay')

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


class Result(DatasetBase):

    def __init__(self, assay):
        self._assay = assay
        self._title = "Clustering analysis of " + assay.sample.title
        self._clusters = None
        self._desc = ""

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

    def _get_clusters(self):
        sample = self.assay.sample

        if not sample.source.has_cell_types:
            return []

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
        self._desc = ""
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


class ProformaWriter(object):

    def __init__(self, builder):
        self._builder = builder
        self._generator = None

    def write(self, project):
        self._write_refeference_proforma(project.reference)
        self._write_common_header(project)
        self._write_field('LC6g', project.title)
        self._write_field('LC6a', project.description)
        self._write_dataset_type(project)
        self._write_cv_terms(project)
        self._write_species(project)
        self._write_count(project)
        self._write_protocols(project)
        self._write_metadata(project)
        self._write_separator()

        for subproject in project.subprojects:
            self._write_subproject(subproject)

        for sample in project.samples:
            self._write_sample(sample)

        self._write_terminator()

    def _write_field(self, field, value):
        self._generator.write_field(field, value)

    def _write_separator(self):
        self._generator.write_separator()

    def _write_terminator(self):
        self._generator.write_terminator()

    def _write_refeference_proforma(self, reference):
        g = self._builder.get_generator(proforma_type=ProformaType.PUB_MINI)
        g.write_header()
        g.write_field('P22', reference.fbrf)
        g.write_field('P2', reference.journal)
        g.write_field('P19')
        g.write_separator()
        self._generator = self._builder.get_generator(proforma_type=ProformaType.DATASET)

    def _write_common_header(self, dataset):
        self._generator.write_header()
        self._write_field('LC1f', 'new')
        self._write_field('LC1a', dataset.symbol)

    def _write_dataset_type(self, dataset):
        self._write_field('LC2a', dataset.entity_type)
        self._write_field('LC2b', dataset.data_type)

    def _write_cv_terms(self, dataset):
        self._write_field('LC13a', dataset.go_cellular_components)
        self._write_field('LC13b', dataset.go_molecular_functions)
        self._write_field('LC13c', dataset.go_biological_processes)
        self._write_field('LC13d', dataset.so_terms)

    def _write_species(self, dataset):
        self._write_field('LC4a', dataset.species)
        self._write_field('LC4i', dataset.other_species)

    def _write_protocols(self, dataset):
        self._write_field('LC11m', dataset.fbcv)
        self._write_field('LC11a', dataset.collection_protocol)
        self._write_field('LC6b', dataset.preparation_protocol)
        self._write_field('LC11c', dataset.assay_protocol)
        self._write_field('LC11e', dataset.analysis_protocol)

    def _write_count(self, dataset):
        self._write_field('LC6d', 'N')
        if dataset.has_count:
            self._write_field('LC6e', dataset.count)
            self._write_field('LC6f', dataset.count_label)

    def _write_metadata(self, project):
        self._write_field('LC7a', "The EMBL-EBI's Single Cell Expression Atlas provides cell-level annotations, clustering data, raw and normalised read counts, and putative marker genes.")
        for source in project.sources:
            self._write_field('LC99a', source.accession)
            self._write_field('LC99b', 'EMBL-EBI Single Cell Expression Atlas Datasets')
        self._write_field('LC8c', f'[{project.lab.name}]({project.lab.url})')

    def _write_subproject(self, subproject):
        self._write_common_header(subproject)
        self._write_field('LC6g', subproject.title)
        self._write_dataset_type(subproject)
        self._write_field('LC3', subproject.project.symbol)
        self._write_cv_terms(subproject)
        self._write_species(subproject)
        self._write_count(subproject)
        self._write_protocols(subproject)
        self._write_separator()

        for sub in subproject.subprojects:
            self._write_subproject(sub)

        for sample in subproject.samples:
            self._write_sample(sample)

    def _write_sample(self, sample):
        self._write_common_header(sample)
        self._write_field('LC6g', sample.title[0].upper() + sample.title[1:])
        self._write_dataset_type(sample)
        self._write_field('LC3', sample.project.symbol)
        self._write_species(sample)
        self._write_field('LC4h', sample.strain)
        self._write_field('LC4f', sample.genotype)
        self._write_tap(sample.developmental_stage, sample.sex, sample.anatomical_part)
        self._write_count(sample)
        self._write_field('LC11m', sample.fbcv)
        self._write_field('LC11a', sample.collection_protocol)
        self._write_separator()
        self._write_assay(sample.assay)

    def _write_assay(self, assay):
        self._write_common_header(assay)
        self._write_field('LC6g', assay.title)
        self._write_dataset_type(assay)
        self._write_field('LC3', assay.sample.project.symbol)
        self._write_field('LC14a', assay.sample.symbol)
        self._write_species(assay.sample)
        self._write_count(assay)
        self._write_field('LC11m', assay.fbcv)
        self._write_separator()
        self._write_result(assay.result)

    def _write_result(self, result):
        self._write_common_header(result)
        self._write_field('LC6g', result.title)
        self._write_dataset_type(result)
        self._write_field('LC3', result.assay.sample.project.symbol)
        self._write_field('LC14b', result.assay.symbol)
        self._write_field('LC14h', 'Dmel R6.32')  # FIXME: Do not hardcode
        self._write_species(result.assay.sample)
        self._write_count(result)
        self._write_separator()

        for cluster in result.clusters:
            self._write_cluster(cluster)

    def _write_cluster(self, cluster):
        self._write_common_header(cluster)
        self._write_field('LC6g', cluster.title)
        self._write_dataset_type(cluster)
        self._write_field('LC3', cluster.result.symbol)
        self._write_species(cluster.result.assay.sample)
        self._write_tap(cluster.result.assay.sample.developmental_stage,
                        cluster.result.assay.sample.sex,
                        cluster.cell_type)
        self._write_count(cluster)
        self._write_separator()

    def _write_tap(self, fbcv, sex, fbbt):
        if sex is not None:
            fbcv += ' | ' + sex
        self._write_field('LC4g', f'<e><t>{fbcv}<a>{fbbt}<s><note>')


class CurationContext(object):

    def __init__(self, config, store, no_exclude, min_cluster_size, with_reads):
        self._store = store
        self._no_exclude = no_exclude
        self._min_cluster = min_cluster_size
        self._with_reads = with_reads
        self._proformae_dir = config.get('curation', 'proformae')

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

    def dataset_from_specfile(self, specfile, with_reads=True):
        spec = json.load(specfile)
        self._expand_defaults(spec)
        dataset = Project(spec, self._min_cluster)

        for source in spec['sources']:
            s = self.source_dataset_from_specfile(source)
            s.apply_corrections()
            dataset.sources.append(s)

        if self._with_reads and with_reads:
            dataset.extract_reads()

        return dataset

    def _expand_defaults(self, spec, parent={}):
        defaults = spec.get('defaults', {})
        self._copy_dict(parent, defaults)

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
@click.option('--with-reads/--without-reads', default=True,
              help="Extract number of reads per biosample.")
@click.pass_context
def curate(ctx, no_exclude, min_cluster_size, with_reads):
    """Access the curation commands."""

    ctx.obj = CurationContext(ctx.obj.config, ctx.obj.raw_store, no_exclude,
                              min_cluster_size, with_reads)
    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-curate> ")
        shell.cmdloop()


@curate.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def extract(ctx, specfile, output):
    """Extract curation data from a dataset."""

    ds = ctx.dataset_from_specfile(specfile)

    for sample in ds.get_all_samples():
        output.write(f"Sample {sample.symbol}\n")
        output.write(f"  Cells: {sample.assay.result.count}\n")
        output.write(f"  Reads: {sample.assay.count}\n")
        for cluster in sample.assay.result.clusters:
            output.write(f"    {cluster.cell_type}: {cluster.count}\n")


@curate.command()
@click.argument('specfile', type=click.File('r'))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.option('--header', '-H', is_flag=True, default=False,
              help="Writes an uncommented header line.")
@click.pass_obj
def sumexpr(ctx, specfile, output, header):
    """Summarize expression data from a dataset."""

    ds = ctx.dataset_from_specfile(specfile, with_reads=False)
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
    builder = ctx.get_proforma_builder(output)
    writer = ProformaWriter(builder)
    writer.write(ds)


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

