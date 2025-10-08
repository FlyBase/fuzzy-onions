# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2025 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import logging
from datetime import datetime


class BaseDatasetExporter(object):
    """Base helper class to export a dataset to another format."""
    
    def __init__(self, db, ontologies, fbbt_corrections):
        self._db = db
        fbbt = ontologies.fbbt.backend
        fbdv = ontologies.fbdv.backend
        self._fbbt_terms = {}
        self._fbdv_terms = {}
        for t in fbdv.terms():
            self._fbdv_terms[t.name] = t.id
        for t in fbbt.terms():
            self._fbbt_terms[t.name] = t.id
            for syn in [s for s in t.synonyms if s.scope == 'EXACT']:
                self._fbbt_terms[syn.description] = t.id
        self._fbbt_corrections = {}
        if fbbt_corrections is not None:
            self._load_fbbt_corrections(fbbt_corrections)
        self._id_by_symbol = {}
        
    def export_as_json(self, datasets):
        """Export all the given datasets into a JSON object."""
        pass
        
    def _get_dataset_id(self, symbol):
        """Get the FBlc identifier for a given dataset symbol."""
        
        did = self._id_by_symbol.get(symbol)
        if did is not None:
            return did
        
        query = f'''SELECT uniquename
                    FROM
                             library
                    WHERE
                             name LIKE '{symbol}';'''
        self._db.cursor.execute(query)
        try:
            did = "FB:" + self._db.cursor.fetchone()[0]
            self._id_by_symbol[symbol] = did
            return did
            return "FB:" + self._db.cursor.fetchone()[0]
        except:
            logging.warning(f"Unknown symbol: {symbol}")
            return "unknown"
        
    def _get_pmid(self, fbrf):
        """Get the PMID for a given FBrf."""

        query = f'''SELECT dbxref.accession
                    FROM
                             pub
                        JOIN pub_dbxref USING (pub_id)
                        JOIN dbxref     USING (dbxref_id)
                        JOIN db         USING (db_id)
                    WHERE
                             db.name LIKE 'pubmed'
                         AND pub.uniquename LIKE '{fbrf}';'''
        self._db.cursor.execute(query)
        try:
            return "PMID:" + self._db.cursor.fetchone()[0]
        except:
            logging.warning(f"Unknown FBrf: {fbrf}")
            return "unknown"
        
    def _propagate_protocols(self, dataset):
        if dataset.collection_protocol:
            for sample in dataset.get_all_samples():
                if not sample.collection_protocol:
                    sample._prot_collection = dataset.collection_protocol
        if dataset.analysis_protocol:
            for result in dataset.get_all_results():
                if not result.analysis_protocol:
                    result._prot_analysis = dataset.analysis_protocol
        
    def _export_fbcv(self, fbcv_value):
        if fbcv_value is None:
            return None
        elif isinstance(fbcv_value, list):
            return [i.split(';')[1].strip() for i in fbcv_value]
        else:
            return fbcv_value.split(';')[1].strip()

    def _export_fbdv(self, term):
        try:
            return self._fbdv_terms[term]
        except KeyError:
            logging.warning(f"Unknown FBdv term: {term}")
            return "FBdv:UNKNOWN"

    def _export_fbbt(self, term):
        term = self._fbbt_corrections.get(term, term)
        try:
            return self._fbbt_terms[term]
        except KeyError:
            logging.warning(f"Unknown FBbt term: {term}")
            return "FBbt:UNKNOWN"

    def _load_fbbt_corrections(self, correction_file):
        with open(correction_file, 'r') as f:
            for line in f:
                old, new = line.strip().split('\t')
                self._fbbt_corrections[old] = new


class BackupDatasetExporter(BaseDatasetExporter):
    """A helper class to export a simplified view of a scRNAseq dataset."""

    def __init__(self, db, ontologies, fbbt_corrections):
        BaseDatasetExporter.__init__(self, db, ontologies, fbbt_corrections)
        
    def export_as_json(self, datasets):
        return {'@type': 'DatasetGroup', 'datasets': [self._export(d) for d in datasets]}

    def _export(self, dataset):
        self._propagate_protocols(dataset)
        d = self._get_new_dict(dataset)
        d['reference'] = self._get_pmid(dataset.reference.fbrf)
        d['study_cvterms'] = [self._export_fbcv(t) for t in dataset.fbcv]
        if dataset.lab:
            d['creator'] = {'name': dataset.lab.name, 'url': dataset.lab.url}
        d['accessions'] = [s.full_accession for s in dataset.sources]
        d['samples'] = self._get_samples(dataset)
        analysis = self._get_project_analysis(dataset)
        if analysis is not None:
            d['analysis'] = analysis
        return d
        
    def _get_new_dict(self, obj):
        return {
            'symbol': obj.symbol,
            'id': self._get_dataset_id(obj.symbol),
            'title': obj.title
            }

    def _add(self, dictionary, name, value):
        if value:
            dictionary[name] = value

    def _get_samples(self, project):
        samples = []
        if project.has_samples:
            results = project.results
            for sample in project.samples:
                s = self._export_sample(sample)
                result = None
                for r in results:
                    if len(r.assays) == 1 and r.assay.sample.symbol == sample.symbol:
                        result = r
                if result is not None and sample.include_clusters:
                    s['analysis'] = self._export_analysis(result)
                samples.append(s)
        elif project.has_subprojects:
            for subproject in project.subprojects:
                sp = self._export_subproject(subproject)
                samples.append(sp)
        return samples

    def _get_project_analysis(self, project):
        extra_results = project.get_extra_results()
        if len(extra_results) == 1:
            a = self._export_analysis(extra_results[0])
            return a
        return None

    def _export_subproject(self, subproject):
        d = self._get_new_dict(subproject)
        d['samples'] = self._get_samples(subproject)
        analysis = self._get_project_analysis(subproject)
        self._add(d, 'analysis', analysis)
        return d

    def _export_sample(self, sample):
        d = self._get_new_dict(sample)
        d['sample_type'] = self._export_fbcv(sample.data_type)
        d['assay_type'] = self._export_fbcv(sample.assay.data_type)
        self._add(d, 'strain', sample.strain)
        self._add(d, 'genotype', sample.genotype)
        d['tissue'] = self._export_fbbt(sample.anatomical_part)
        d['stage'] = self._export_fbdv(sample.developmental_stage)
        self._add(d, 'sex', sample.sex)
        if sample.assay.technical_reference:
            d['technical_reference'] = self._get_dataset_id(
                sample.assay.technical_reference
            )
        if sample.assay.biological_reference:
            d['biological_reference'] = self._get_dataset_id(
                sample.assay.biological_reference
            )
        if sample.entities:
            entities = []
            for entity, entity_type in sample.entities:
                entities.append({'entity': "FB:" + entity, 'entity_type': entity_type})
            d['experimental_factors'] = entities
        self._add(d, 'collection_cvterms', self._export_fbcv(sample.fbcv))
        self._add(d, 'assay_cvterms', self._export_fbcv(sample.assay.fbcv))
        self._add(d, 'collection_protocol', sample.collection_protocol)
        return d

    def _export_analysis(self, result):
        d = self._get_new_dict(result)
        d['analysis_type'] = self._export_fbcv(result.data_type)
        self._add(d, 'analysis_protocol', result.analysis_protocol)
        self._add(d, 'analysis_cvterms', self._export_fbcv(result.fbcv))
        d['cell_count'] = result.count
        clusters = []
        for cluster in result.clusters:
            c = self._get_new_dict(cluster)
            c['cell_type'] = self._export_fbbt(cluster.cell_type)
            c['cell_count'] = cluster.count
            clusters.append(c)
        d['clusters'] = clusters
        return d
    

class AllianceDatasetExporter(BaseDatasetExporter):
    """A helper class to export datasets in the format expected by the Alliance."""
    
    def __init__(self, db, ontologies, fbbt_corrections):
        BaseDatasetExporter.__init__(self, db, ontologies, fbbt_corrections)
        self._today = datetime.now().isoformat(timespec='seconds')
        
    def export_as_json(self, datasets):
        dataset_objects = []
        sample_objects = []
        for dataset in datasets:
            self._export_project(dataset, dataset_objects, sample_objects)
        return {
            "high_throughput_expression_dataset_annotation_ingest_set": dataset_objects,
            "high_throughput_expression_dataset_sample_annotation_ingest_set": sample_objects
            }
    
    def _get_new_dict(self):
        # FIXME: Should we add a "created_by_curie"? If so how do we get the
        # curie of the curator?
        return {
            "date_created": self._today,
            "internal": False,
            "obsolete": False
            }
        
    def _export_project(self, project, dataset_objects, sample_objects):
        if project.is_top_project:
            self._propagate_protocols(project)

        d = self._get_new_dict()
        d["primary_external_id"] = self._get_dataset_id(project.symbol)
        # FIXME: The dataset ID in "htp_expression_dataset_dto" is the same as
        # the primary_external_id (and is in fact the MOD ID), because apart
        # from the top-level project, none of those objects correspond to an
        # external database entity.
        d["htp_expression_dataset_dto"] = {"curie": self._get_dataset_id(project.symbol)}
        d["name"] = project.title
        # FIXME: The "symbol" slot is currently only defined on the
        # HTPExpressionDatasetAnnotation, _not_ on its *DTO counterpart, so
        # the following is not compliant.
        d["symbol"] = project.symbol
        
        if project.is_top_project:
            d["htp_expression_dataset_dto"]["secondary_identifiers"] = [s.full_accession for s in project.sources]
            d["references_curies"] = ["FB:" + project.reference.fbrf]
            # FIXME: Translate from FBcv terms to "Data Set Category Tags"
            d["category_tag_names"] = ["TODO"]
            
        if project.has_subprojects:
            # FIXME: Check that it is OK to use sub_series for this
            d["sub_series"] = []
            for subproject in project.subprojects:
                d["sub_series"].append(self._get_dataset_id(subproject.symbol))
                self._export_project(subproject, dataset_objects, sample_objects)
                
        for sample in project.samples:
            sample_objects.append(self._export_sample(sample, project))
            
        extra_results = project.get_extra_results()
        if len(extra_results) == 1:
            # FIXME: For now it is not possible to have an analysis object at
            # the level of the dataset (only at the level of samples), so the
            # following is not compliant
            d["analysis"] = self._export_analysis(extra_results[0])
            
        dataset_objects.append(d)
            
        
    def _export_sample(self, sample, parent):
        d = self._get_new_dict()
        d["primary_external_id"] = self._get_dataset_id(sample.symbol)
        d["htp_expression_sample_dto"] = {
            # FIXME: Same issue as for project objects (duplicate IDs)
            "curie": self._get_dataset_id(sample.symbol)
            }
        d["htp_expression_sample_title"] = sample.title
        # FIXME: The "symbol" slot is missing from the
        # HTPExpressionDatasetSampleAnnotationDTO class, so the following is not
        # compliant. 
        d["symbol"] = sample.symbol
        d["dataset_ids"] = [self._get_dataset_id(parent.symbol)]
        # FIXME: Shouldn't that slot be "taxon_curie" in the DTO class?
        d["taxon"] = "NCBITaxon:7227"
        # FIXME: Shouldn't that slot be "genetic_sex_name" in the DTO class?
        d["genetic_sex"] = sample.sex or "pooled sexes"
        # FIXME: Translate from Fbcv to OBI
        d["htp_expression_sample_type_curie"] = "TODO"
        # FIXME: Translate from FBcv to MMO
        d["expression_assay_used"] = "TODO"
        d["htp_expression_sample_location_dtos"] = [{
            "anatomical_structure_curie": self._export_fbbt(sample.anatomical_part)
            }]
        d["htp_expression_sample_age_dto"] = {
            "stage_dto": {
                "developmental_stage_start_curie": self._export_fbdv(sample.developmental_stage)
                }
            }
        if sample.genotype is not None:
            # FIXME: Should ideally use AGM models instead (much more complex
            # but would allow to represent the exact genetic entities involved,
            # when they are known).
            d["genomic_information_dto"] = {
                "biosample_text": sample.genotype
                }
            
        
        result = None
        for r in parent.results:
            if len(r.assays) == 1 and r.assay.sample.symbol == sample.symbol:
                result = r
        if result is not None and sample.include_clusters:
            # FIXME: There is currently no "analysis" slot in the
            # HTPExpressionDatasetSampleAnnotation DTO class (there is only a
            # "analysis" slot in the non-DTO counterpart class), so the
            # following is not compliant.
            d["analysis"] = self._export_analysis(result)
            
        return d
            
    
    def _export_analysis(self, result):
        # FIXME: Shouldn't there be a AnalysisDTO class?
        # With something like a "analysis_type_curie" slot?
        d = self._get_new_dict()
        d["primary_external_id"] = self._get_dataset_id(result.symbol)
        d["analysis_type"] = self._export_fbcv(result.data_type)
        d["cell_count"] = result.count
        d["symbol"] = result.symbol
        if result.analysis_protocol is not None:
            d["note_dtos"] = [{
                "free_text": result.analysis_protocol,
                 "note_type_name": "analysis_protocol"
                }]
        clusters = []
        for cluster in result.clusters:
            # FIXME: Shouldn't there be a SingleCellExpressionClusterDTO class?
            # With something like a "cell_type_curie" slot?
            c = self._get_new_dict()
            c["primary_external_id"] = self._get_dataset_id(cluster.symbol)
            c["symbol"] = cluster.symbol
            c["title"] = cluster.title
            c["cell_count"] = cluster.count
            c["cell_type"] = self._export_fbbt(cluster.cell_type)
            clusters.append(c)
        d["cell_clusters"] = clusters
        return d
