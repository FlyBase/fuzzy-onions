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
import yaml
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
    """A FlyBase paper."""

    def __init__(self, spec):
        self._fbrf = spec.get('fbrf', None)
        self._journal = spec.get('journal', None)

    @property
    def fbrf(self):
        """The FlyBase identifier for the paper."""
        
        return self._fbrf

    @property
    def journal(self):
        """The abbreviated form of the journal."""
        
        # TODO: Fetch from the database with the FBrf?
        return self._journal


class SourceLab(object):
    """The research group or organisation that created the dataset.
    
    This feeds the LC8c field of the dataset proforma.
    """

    def __init__(self, spec):
        self._name = spec.get('name', None)
        self._url = spec.get('url', None)

    @property
    def name(self):
        """Human-readable name of the research group."""
        
        return self._name

    @property
    def url(self):
        """URL to the public web site of the research group."""
        
        return self._url


class DatasetBase(object):
    """Base class for all dataset objects.
    
    This class holds data fields that are common to most dataset
    objects, regardless of the level.
    """

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
        """The dataset symbol, directly usable in LC1a."""
                
        return self._symbol

    @property
    def title(self):
        """The dataset title, directly usable in LC6g."""
        
        return self._title

    @property
    def description(self):
        """A long description of the dataset (LC6a, optional)."""
        
        if hasattr(self, '_desc'):
            return self._desc
        else:
            return None

    @property
    def entity_type(self):
        """The type of dataset object (LC2a).
        
        This property is overriden in subclasses to automatically
        return the correct value.
        """
        
        pass

    @property
    def data_type(self):
        """The type of data represented in the dataset object (LC2b).
        
        This property is overriden in subclasses to automatically
        return the correct value.
        """
        
        pass

    @property
    def go_cellular_components(self):
        """The cellular components GO terms (LC13a).
        
        This property is an array of strings already in the format
        expected in the LC13a field ("term name ; GO:xxxxxxx"). It
        may be None if GO curation has not been performed, or an
        empty array if the curators explicitly assigned no GO terms.
        """
        
        return self._go_cc

    @property
    def go_molecular_functions(self):
        """The molecular functions GO terms (LC13b).
        
        This property is an array of strings already in the format
        expected in the LC13b field ("term name ; GO:xxxxxxx"). It
        may be None if GO curation has not been performed, or an
        empty array if the curators explicitly assigned no GO terms.
        """
        
        return self._go_mf

    @property
    def go_biological_processes(self):
        """The biological processes GO terms (LC13c).
        
        This property is an array of strings already in the format
        expected in the LC13c field ("term name ; GO:xxxxxxx"). It
        may be None if GO curation has not been performed, or an
        empty array if the curators explicitly assigned no GO terms.
        """
        
        return self._go_bp

    @property
    def so_terms(self):
        """The Sequence Ontology terms (LC13d).
        
        This property is an array of strings already in the format
        expected in the LC13d field ("term name ; SO:xxxxxxx"). It
        may be None if SO curation has not been performed, or an
        empty array if the curators explicitly assigned no SO terms.
        """
        return self._so

    @property
    def fbcv(self):
        """The FlyBase Controlled Vocabulary terms (LC11m).
        
        This property is an array of strings already in the format
        expected in the LC11m field ("term name ; FBcv:xxxxxxx"). It
        may be None if FBcv curation has not been performed, or an
        empty array if the curators explicitly assigned no FBcv
        terms.
        """
        
        return self._fbcv

    @property
    def collection_protocol(self):
        """The sample collection or isolation protocol (LC11a)."""
        
        return self._prot_collection

    @property
    def preparation_protocol(self):
        """The sample preparation protocol (LC6b)."""
        
        return self._prot_preparation

    @property
    def assay_protocol(self):
        """The protocol for the assay proper (LC11c)."""
        
        return self._prot_assay

    @property
    def analysis_protocol(self):
        """The protocol for data analysis (LC11e)."""
        
        return self._prot_analysis

    @property
    def species(self):
        """The species of derivation (LC4a).
        
        This property contains a 4-letter species code (typically
        Dmel).
        """
        
        return self._species

    @property
    def other_species(self):
        """A list of extra derivation species (LC4i)."""
        
        return []

    @property
    def has_count(self):
        """Indicate whether the dataset has a numerical component.
        
        If set to True for a given dataset object, the `count`
        property of that object will contain the value expected for
        the LC6e field.
        """
        
        return False

    @property
    def count_label(self):
        """Label for the entities count value (LC6f)."""
        
        return None


class ProjectContainer(DatasetBase):
    """A dataset object that can contains subprojects or samples.
    
    An object of this type will either contain subprojects, or it
    will contain samples. It is not expected that a project will
    contain both subprojects and samples.
    """

    def __init__(self, spec, parent=None):
        super().__init__(spec)

        self._parent = parent

        self._subprojects = []
        for subproject in spec.get('subprojects', []):
            self._subprojects.append(ProjectContainer(subproject, self))

        self._samples = []
        for sample in spec.get('samples', []):
            self._samples.append(Biosample(sample, self))

        self._results = []
        for sample in self._samples:
            self._results.append(Result(sample.assay))

        if 'result' in spec:
            symbol = spec['result'].get('symbol')
            title = spec['result'].get('title')
            stage = spec['result'].get('stage')
            r = Result([s.assay for s in self._samples],
                       project=self,
                       symbol=symbol,
                       title=title,
                       stage=stage)
            self._results.append(r)

    @property
    def subprojects(self):
        """The subprojects contained in this project (may be empty).
        
        :return: a list of :class:`ProjectContainer` objects
        """
        
        return self._subprojects

    @property
    def has_subprojects(self):
        """Indicate whether this project contains subprojects."""
        
        return len(self.subprojects) > 0

    @property
    def samples(self):
        """The samples contained in this project (may be empty).
        
        :return: a list of :class:`Biosample` objects
        """
        
        return self._samples

    @property
    def has_samples(self):
        """Indicate whether this project contains samples."""
        
        return len(self.samples) > 0

    @property
    def project(self):
        """The (sub)project this subproject belongs to.
        
        :return: a :class:`ProjectContainer` object, or None if the
            current project is a top-level project
        """
        
        return self._parent

    @property
    def is_top_project(self):
        """Indicate whether the current project is a top-level project."""
        
        return self._parent is None

    @property
    def top_project(self):
        """The top-level project this (sub)project is a part of.
        
        If the current object is already a top-level project, this
        property returns the current object itself.
        
        :return: a :class:`Project` object
        """
        
        p = self
        while not p.is_top_project:
            p = self.project
        return p

    @property
    def results(self):
        """All the results for the current project.
        
        If the current project contains samples, then this property
        returns the results for all the samples in the project, in
        addition to any project-defined result. If it contains
        subprojects, then this only returns the results defined by
        the current project; use the :meth:`get_all_results` method
        to get the results from all the subprojects.
        
        :return: a list of :class:`Result` objects
        """
        
        return self._results

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
        """Get all samples from all subprojects.
        
        :return: a list of :class:`Biosample` objects
        """
        
        if self.has_samples:
            return self.samples
        else:
            samples = []
            for subproject in self.subprojects:
                samples.extend(subproject.get_all_samples())
            return samples

    def get_all_results(self):
        """Get all results from all subprojects.
        
        :return: a list of :class:`Result` objects
        """
        
        results = []
        results.extend(self.results)
        for subproject in self.subprojects:
            results.extend(subproject.get_all_results())
        return results

    def get_assay_single_analysis(self, assay):
        """Get the assay-specific result for the specified assay.
        
        :return: a :class:`Result` object
        """
        
        results = [r for r in self.results if r.is_single_analysis_of(assay)]
        if len(results) == 1:
            return results[0]
        else:
            # This should not happen
            return None

    def get_extra_results(self):
        """Get all project-defined results.
        
        This method returns all the results that are *not*
        assay-specific results.
        
        :return: a list of :class:`Result` objects
        """
        
        assays = [s.assay for s in self.samples]
        extra = []
        for result in self.results:
            is_extra = True
            for assay in assays:
                if result.is_single_analysis_of(assay):
                    is_extra = False
            if is_extra:
                extra.append(result)
        return extra


class Project(ProjectContainer):
    """A dataset object representing a top-level project."""

    def __init__(self, spec, min_cluster_size=0):
        super().__init__(spec)
        self._min_cluster_size = min_cluster_size

        self._reference = FlybaseReference(spec.get('reference', {}))
        self._lab = SourceLab(spec.get('creator', {}))
        self._ref_genome = spec.get('reference_genome', 'Dmel R6.32')

        self._sources = []

    @property
    def sources(self):
        """The list of EBI source datasets for this project.
        
        :return: a list of :class:`SourceDataset` objects
        """
        
        return self._sources

    @property
    def other_species(self):
        species = []
        for sample in self.get_all_samples():
            species.extend(sample.other_species)
        return species

    @property
    def reference(self):
        """The FlyBase paper reporting this dataset.
        
        :return: a :class:`FlyBaseReference` object
        """
        
        return self._reference

    @property
    def lab(self):
        """The research group that created this dataset.
        
        :return: a :class:`SourceLab` object
        """
        
        return self._lab

    @property
    def min_cluster_size(self):
        """The minimal number of cells in a cluster.
        
        If a cluster in the source dataset contains less than this
        threshold, the cluster is excluded from the results.
        """
        
        return self._min_cluster_size

    @property
    def reference_genome(self):
        """The version of the reference genome assembly (LC14h)."""
        
        return self._ref_genome

    def extract_reads(self):
        """Extract the number of sequencing reads for each sample.
        
        This method parses the EBI "Raw Expression Matrix" to get the
        number of sequencing reads for each sample described in this
        project (or its subprojects).
        
        After extraction, the number of sequencing reads is available
        in the :attr:`count` property of the :class:`Assay` object
        associated with each sample.
        """
        
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
        """Produce the Summarised Expression Table for the dataset.
        
        This method parses the EBI "Normalised Expression Matrix" to
        compute, for each gene in each cluster for each sample, the
        average expression of that gene in that cluster, and the
        "extent of expression" (the proportion of cells of that
        cluster in which the gene is expressed).
        
        :return: a :class:`pandas.DataFrame` object
        """
        
        # Pre-processing
        # We prepare a data structure for each cluster in each sample.
        # That structure will be filled as we read the expression table
        # below. This allows to read the table only once.

        clusters_by_cell_id = {}
        results = self.get_all_results()
        for result in results:
            for cluster in result.clusters:
                for cell_id in cluster.get_cell_ids():
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

        table = None
        for result in results:
            for cluster in result.clusters:
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
                if table is not None:
                    table = table.append(d)
                else:
                    table = d

        logging.info("Post-processing complete.")

        return table


class Biosample(DatasetBase):
    """A biological sample in a dataset."""

    def __init__(self, spec, parent):
        super().__init__(spec)
        self._project = parent
        if self._desc is None:
            self._desc = ""

        self._is_sn = spec.get('is_single_nucleus', False)

        self._other_species = spec.get('species', {}).get('other', [])

        self._strain = spec.get('strain', None)
        self._genotype = spec.get('genotype', None)

        self._fbbt = spec.get('anatomical_part', None)
        self._fbdv = spec.get('stage', None)
        self._sex = spec.get('sex', None)

        self._fbcv = spec.get('cv_terms', {}).get('fbcv_sample')

        self._entities = spec.get('entities', [])

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
        """The project this sample belongs to (LC3).
        
        :return: a :class:`ProjectContainer` object
        """
        
        return self._project

    @property
    def top_project(self):
        """The top-level this sample belongs to.
        
        If the dataset does not contain any subproject, then this
        property returns the same value as :attr:`project`.
        
        :return: a :class:`Project` object
        """
        
        return self.project.top_project

    @property
    def is_single_nucleus(self):
        """Indicate whether the sample contains single nuclei."""
        
        return self._is_sn

    @property
    def other_species(self):
        return self._other_species

    @property
    def strain(self):
        """The fly strain that provided this sample (LC4h)."""
        
        return self._strain

    @property
    def genotype(self):
        """The genotype of the flies that provided this sample (LC4f)."""
        
        return self._genotype

    @property
    def anatomical_part(self):
        """The organ or tissue the sample comes from.
        
        This is a term from the Drosophila Anatomy Ontology (FBbt).
        """
        
        return self._fbbt

    @property
    def developmental_stage(self):
        """The stage at which the sample has been collected.
        
        This is a term from the Drosophila Developmental Ontology
        (FBdv).
        """
        
        return self._fbdv

    @property
    def sex(self):
        """The sex of the individuals that provided the sample.
        
        Can either be `female`, `male`, or None if either the sex is
        unknown or the sample comes from a mix of both female and
        male individuals.
        """
        
        return self._sex

    @property
    def entities(self):
        """The list of experimental entities involved.
        
        This is a list of tuple, where each tuple contains the
        FlyBase identifier for the experimental entity (LC12a) and
        a value indicating the type of entity (LC12b).
        """
        
        return self._entities

    @property
    def assay(self):
        """The assay performed on the sample.
        
        :return: a :class:`Assay` object
        """
        
        return self._assay

    @property
    def subset(self):
        """The cells that make up the sample.
        
        This property returns a filtered view of the EBI's original
        "Experiment Design Table", containing only the cells that
        belong to the current sample.
        """
        
        if self._subset is None:
            self._subset = self._get_subset()
        return self._subset

    @property
    def source(self):
        """The source datasets with the data for this sample.
        
        :return: a list of :class:`SourceDataset` object
        """
        
        if self._source is None:
            self._source = self._get_source()
        return self._source

    def exclude_cell_type(self, cell_type):
        """Indicate whether a cell type should be excluded.
        
        This method takes a single cell type (as a FBbt term) and
        indicates whether, according to curators' instructions
        (either defined at the level of the current sample or at any
        higher level), that cell type should be excluded from any
        further analysis.
        
        :return: True if the cell type shall be excluded
        """
        
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
    """A single-cell RNA-Seq assay."""

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
        """The assay that is a technical reference to this one (LC14e).
        
        :return: a :class:`Assay` object
        """
        
        return self._get_full_reference(self._tech_ref)

    @property
    def biological_reference(self):
        """The assay that is a biological reference to this one (LC14f).
        
        :return: a :class:`Assay` object
        """
        
        return self._get_full_reference(self._biol_ref)

    @property
    def sample(self):
        """The sample from which the assay is derived (LC14a).
        
        :return: a :class:`Biosample` object
        """
        
        return self._sample

    def _get_full_reference(self, reference):
        if reference is None or len(reference) == 0:
            return reference

        symbol = self.sample.project.top_project.symbol + '_' + reference
        samples = [s for s in self.sample.project.top_project.get_all_samples() if s.symbol == symbol]
        if len(samples) == 1:
            return samples[0].assay.symbol
        elif len(samples) == 0:
            logging.warn(f"Referenced assay not found: {reference}")
        elif len(samples) > 1:
            logging.warn(f"Assay reference ambiguous: {reference}")
        return reference


class Result(DatasetBase):
    """A clustering analysis from a scRNAseq assay."""

    def __init__(self, assay, project=None, symbol=None, title=None, stage=None):
        if isinstance(assay, list):
            self._assays = assay
        else:
            self._assays = [assay]

        if project is None:
            # This is a single assay-specific analysis
            project = self._assays[0].sample.project
            symbol = self._assays[0].symbol + '_clustering'
            title = "Clustering analysis of " + self._assays[0].sample.title
        else:
            if symbol is None:
                symbol = project.symbol + '_clustering'
            else:
                symbol = project.symbol + '_' + symbol
            if title is None:
                title = "Clustering analysis of " + project.title[0].lower() + project.title[1:]

        self._project = project
        self._symbol = symbol
        self._title = title
        self._stage = stage

        self._desc = ""
        self._clusters = None

    @property
    def entity_type(self):
        return 'result ; FBcv:0003026'

    @property
    def data_type(self):
        return 'cell clustering analysis ; FBcv:0009002'

    @property
    def assay(self):
        """The assay from which this result is derived (LC14b).
        
        This property returns the first assay from which the current
        analysis is derived. Note that an analysis may derive from
        several assays. The :attr:`assays` property returns all
        those assays and its use should be preferred.
        
        :return: a :class:`Assay` object
        """
        
        return self._assays[0]

    @property
    def assays(self):
        """The assays from which this result is derived (LC14b).
        
        :return: a list of :class:`Assay` objects
        """
        
        return self._assays

    @property
    def project(self):
        """The project this analysis belongs to (LC3).
        
        :return: a :class:`ProjectContainer` object
        """
        
        return self._project

    @property
    def analysis_protocol(self):
        # FIXME: Allow to specify custom protocol for
        # non-sample-specific analyses
        return self.assay.sample.analysis_protocol

    @property
    def count(self):
        n = 0
        for assay in self.assays:
            n += len(assay.sample.subset)
        return n

    @property
    def has_count(self):
        return True

    @property
    def count_label(self):
        return 'cells'

    @property
    def clusters(self):
        """The clusters identified in this analysis.
        
        :return: a list of :class:`Cluster` objects
        """
        
        if self._clusters is None:
            self._clusters = self._get_clusters()
        return self._clusters

    def get_sex(self):
        """Get the sex associated with this analysis.
        
        This method returns the sex common to all samples from which
        this analysis is derived, if any. If the analysis is derived
        from samples of different sexes, None is returned.
        
        :return: `female`, `male`, or None
        """
        
        all_sexes = set([a.sample.sex for a in self.assays])
        if len(all_sexes) == 1:
            return all_sexes.pop()
        else:
            return None

    def is_single_analysis_of(self, assay):
        """Indicate whether this analysis is derived from the specified assay."""
        
        if len(self.assays) > 1:
            return False
        return self.assay.symbol == assay.symbol

    def _get_clusters(self):
        clusters = []
        cell_types = {}
        sex = self.get_sex()

        for assay in self.assays:
            sample = assay.sample
            if not sample.source.has_cell_types:
                continue

            for cell_type in sample.subset[sample.source.cell_type_column].dropna().unique():
                if sample.exclude_cell_type(cell_type):
                    continue

                n = len(sample.subset.loc[sample.subset[sample.source.cell_type_column] == cell_type])
                if cell_type in cell_types:
                    cell_types[cell_type][0] += n
                else:
                    cell_types[cell_type] = [n, sample]

        threshold = self.project.top_project.min_cluster_size
        for cell_type in cell_types.keys():
            n, sample = cell_types[cell_type]
            if n >= threshold:
                clusters.append(Cluster(cell_type, n, self, sample, self._stage, sex))

        return clusters


class Cluster(DatasetBase):
    """A transcriptional cell cluster."""

    def __init__(self, cell_type, count, result, sample, stage, sex):
        self._cell_type = cell_type
        self._count = count
        self._result = result
        self._sample = sample
        self._simple_ct = None
        self._stage = stage
        self._sex = sex
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
        """The clustering analysis that identified this cluster.
        
        :return: a :class:`Result` object
        """
        
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
        """The type of cells in this cluster.
        
        This is a FBbt term.
        """
        
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
        """A simplified version of the cell type.
        
        The simplified cell type is a version of the cell type that
        is intended to be used in places where there is limited space
        and full cell type names could be too long.
        
        The simplified cell type is either specified directly by the
        curators, or derived from the application of a few
        simplifying rules.
        """
        
        if self._simple_ct is None:
            sct = self._sample.source.get_simplified_cell_type(self._cell_type)
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
    def developmental_stage(self):
        if self._stage is None:
            return self._sample.developmental_stage
        else:
            return self._stage

    @property
    def sex(self):
        return self._sex

    def get_cell_ids(self):
        """Get the IDs of all cells in this cluster."""
        
        ids = []

        for assay in self.result.assays:
            sample = assay.sample
            if not sample.source.has_cell_types:
                continue

            ct_column = sample.source.cell_type_column
            cl_subset = sample.subset.loc[sample.subset[ct_column] == self.cell_type]
            ids.extend(cl_subset['Assay'].values)

        return ids


class ProformaWriter(object):
    """Helper class to produce a proforma from a dataset."""

    def __init__(self, builder):
        """Create a new instance.
        
        :param builder: a :class:`ProformaGeneratorBuilder` object
        """
        
        self._builder = builder
        self._generator = None

    def write(self, project):
        """Generate a proforma for the specified project.
        
        :param project: a :class:`Project` object
        """
        
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

        for result in project.get_extra_results():
            self._write_result(result)

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

        for result in subproject.get_extra_results():
            self._write_result(result)

    def _write_sample(self, sample):
        self._write_common_header(sample)
        self._write_field('LC6g', sample.title[0].upper() + sample.title[1:])
        self._write_dataset_type(sample)
        self._write_field('LC3', sample.project.symbol)
        self._write_species(sample)
        self._write_field('LC4h', sample.strain)
        self._write_field('LC4f', sample.genotype)
        self._write_tap(sample.developmental_stage, sample.sex, sample.anatomical_part)
        self._write_entities(sample.entities)
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
        self._write_field('LC14d', assay.technical_reference)
        self._write_field('LC14e', assay.biological_reference)
        self._write_species(assay.sample)
        self._write_count(assay)
        self._write_field('LC11m', assay.fbcv)
        self._write_separator()
        self._write_result(assay.sample.project.get_assay_single_analysis(assay))

    def _write_result(self, result):
        self._write_common_header(result)
        self._write_field('LC6g', result.title)
        self._write_dataset_type(result)
        self._write_field('LC3', result.project.symbol)
        self._write_field('LC14b', '\n'.join([a.symbol for a in result.assays]))
        self._write_field('LC14h', result.project.top_project.reference_genome)
        self._write_species(result.assay.sample)
        self._write_count(result)
        self._write_field('LC11e', result.analysis_protocol)
        self._write_separator()

        for cluster in result.clusters:
            self._write_cluster(cluster)

    def _write_cluster(self, cluster):
        self._write_common_header(cluster)
        self._write_field('LC6g', cluster.title)
        self._write_dataset_type(cluster)
        self._write_field('LC3', cluster.result.symbol)
        self._write_species(cluster.result.assay.sample)
        self._write_tap(cluster.developmental_stage, cluster.sex, cluster.cell_type)
        self._write_count(cluster)
        self._write_separator()

    def _write_tap(self, fbcv, sex, fbbt):
        if sex is not None:
            fbcv += ' | ' + sex
        self._write_field('LC4g', f'<e><t>{fbcv}<a>{fbbt}<s><note>')

    def _write_entities(self, entities):
        for entity_ref, entity_type in entities:
            self._write_field('LC12a', entity_ref)
            self._write_field('LC12b', entity_type)


class CurationContext(object):
    """A helper object for all curation operations."""

    def __init__(self, config, store, no_exclude, min_cluster_size, with_reads):
        """Create a new instance.
        
        :param config: the main configuration object
        :param store: the dataset local file store
        :param no_exclude: do not exclude any cell type if True
        :param min_cluster_size: exclude any cluster containing less
            than the specified number of cells
        :param with_reads: if True, extract the number of sequencing
            reads and write the count in the Assay proforma
        """
        
        self._store = store
        self._no_exclude = no_exclude
        self._min_cluster = min_cluster_size
        self._with_reads = with_reads
        self._proformae_dir = config.get('curation', 'proformae')

    def get_proforma_builder(self, output):
        """Get a proforma writer set up to write to the specified file.
        
        :param output: a file-like object
        :return: a :class:`ProformaGeneratorBuilder` object
        """
        
        return ProformaGeneratorBuilder(self._proformae_dir, output)

    def source_dataset_from_specfile(self, specfile):
        """Build a source dataset object from a JSON specification.
        
        :param specfile: the name of a file containing a JSON or YAML
            description of the source dataset
        :return: a :class:`SourceDataset` object
        """

        spec = self.load_spec_file(specfile)
        dataset = SourceDataset(spec)
        dataset.feed_data(self._store)
        return dataset

    def dataset_from_specfile(self, specfile, with_reads=True):
        """Build a FlyBase dataset object from a JSON specification.
        
        :param specfile: the name of a file containing a JSON or YAML
            description of the dataset
        :param with_reads: if True (the default), automatically
            extract the number of sequencing reads when building
            the dataset
        :return: a :class:`Project` object
        """

        spec = self.load_spec_file(specfile)
        self.expand_defaults(spec)
        dataset = Project(spec, self._min_cluster)

        for source in spec['sources']:
            s = self.source_dataset_from_specfile(source)
            s.apply_corrections()
            dataset.sources.append(s)

        if self._with_reads and with_reads:
            dataset.extract_reads()

        return dataset

    def load_spec_file(self, filename):
        """Load a JSON or YAML specification file.
        
        :param filename: the name of the specification file
        :return: a dictionary representing the specification
        """

        is_json = filename.endswith('.json')
        with open(filename, 'r') as f:
            if is_json:
                spec = json.load(f)
            else:
                spec = yaml.load(f, Loader=yaml.SafeLoader)

        return spec

    def expand_defaults(self, spec, parent={}):
        """Propagate default values in a JSON description.
        
        If the JSON description of a project contains a dictionary
        named `defaults`, the values of that dictionary are
        propagated to all samples in the project (and subprojects).
        
        :param spec: the JSON dictionary describing the project
        :param parent: supplementary default values to apply, in
            addition to those that may be defined in the project
        """
        
        defaults = spec.pop('defaults', {})
        self._copy_dict(parent, defaults)

        for subproject in spec.get('subprojects', []):
            self.expand_defaults(subproject, defaults)

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
@click.argument('specfile', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def extract(ctx, specfile, output):
    """Extract curation data from a dataset.
    
    This command is mostly intended for debugging. It reads a dataset
    description file and lists all the samples from the project,
    along with the clusters found in each sample (with the number of
    cells in each sample and each cluster).
    """

    ds = ctx.dataset_from_specfile(specfile)

    for sample in ds.get_all_samples():
        output.write(f"Sample {sample.symbol}\n")

        result = sample.project.get_assay_single_analysis(sample.assay)
        output.write(f"  Cells: {result.count}\n")
        output.write(f"  Reads: {sample.assay.count}\n")
        for cluster in result.clusters:
            output.write(f"    {cluster.cell_type}: {cluster.count}\n")


@curate.command()
@click.argument('specfile', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.option('--header', '-H', is_flag=True, default=False,
              help="Writes an uncommented header line.")
@click.pass_obj
def sumexpr(ctx, specfile, output, header):
    """Summarize expression data from a dataset.
    
    This command produces the Summarised Expression Table from a
    dataset description file. The table is written in the format
    expected by Harvard developers, ready to be uploaded to the
    Harvard FTP server.
    """

    ds = ctx.dataset_from_specfile(specfile, with_reads=False)
    result = ds.summarise_expression()
    if not header:
        # Write a commented header line (needed for harvdev processing)
        output.write('#genes\t')
        output.write('\t'.join(result.columns))
        output.write('\n')
    result.to_csv(output, sep='\t', header=header, float_format='%.6f')


@curate.command()
@click.argument('specfile', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def proforma(ctx, specfile, output):
    """Generate a proforma for a dataset.
    
    This command produces a dataset proforma from a dataset
    description file. Depending on the extent of the description,
    the generated proforma may be directly loadable or may require
    further edits by the curators.
    
    In any case, the produced proforma should always be checked by
    Peeves before being submitted for loading.
    """

    ds = ctx.dataset_from_specfile(specfile)
    builder = ctx.get_proforma_builder(output)
    writer = ProformaWriter(builder)
    writer.write(ds)


@curate.command()
@click.argument('specfile', type=click.Path(exists=True))
@click.pass_obj
def fixscea(ctx, specfile):
    """Generate correction files for the SCEA.
    
    This command takes a source dataset description file and produces
    the correction files needed by the EBI curators, if such
    corrections have been described as necessary.
    """

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


@curate.command()
@click.argument('specfile', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.File('w'), default='-',
              help="Write to the specified file instead of standard output.")
@click.pass_obj
def expand(ctx, specfile, output):
    """Expand a dataset description file.
    
    This command is mostly intended for debugging. It takes a JSON
    file describing a dataset, expands any default values in it, and
    writes a new JSON file with all default values expanded. This
    allows to check that default values are expanded as intended by
    the curators.
    """

    spec = ctx.load_spec_file(specfile)
    ctx.expand_defaults(spec)
    json.dump(spec, output, indent=4)

