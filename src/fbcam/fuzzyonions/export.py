# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2025 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.


class DictDatasetExporter(object):
    """A helper class to export a simplified view of a scRNAseq dataset."""

    def __init__(self, db):
        self._db = db

    def _get_dataset_id(self, symbol):
        """Get the FBlc identifier for a given dataset symbol."""

        query = f'''SELECT uniquename
                    FROM
                             library
                    WHERE
                             name LIKE '{symbol}';'''
        self._db.cursor.execute(query)
        try:
            return self._db.cursor.fetchone()[0]
        except:
            print(f"Unknown symbol: {symbol}")
            return "unknown"

    def export(self, dataset):
        self._propagate_protocols(dataset)
        d = {}
        self._fill_id_slots(dataset, d)
        d['reference'] = dataset.reference.fbrf
        d['study_type'] = [t for t in dataset.fbcv]
        if dataset.lab:
            d['creator'] = {'name': dataset.lab.name, 'url': dataset.lab.url}
        d['accessions'] = [s.accession for s in dataset.sources]
        d['samples'] = self._get_samples(dataset)
        analysis = self._get_project_analysis(dataset)
        if analysis is not None:
            d['analysis'] = analysis
        return d

    def _propagate_protocols(self, dataset):
        if dataset.collection_protocol:
            for sample in dataset.get_all_samples():
                if not sample.collection_protocol:
                    sample._protocol = dataset.collection_protocol
        if dataset.analysis_protocol:
            for result in dataset.get_all_results():
                if not result.analysis_protocol:
                    result._protocol = dataset.analysis_protocol

    def _fill_id_slots(self, obj, dictionary):
        dictionary['symbol'] = obj.symbol
        # dictionary['id'] = self._get_dataset_id(obj.symbol)
        dictionary['title'] = obj.title

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
        d = {}
        self._fill_id_slots(subproject, d)
        d['samples'] = self._get_samples(subproject)
        analysis = self._get_project_analysis(subproject)
        self._add(d, 'analysis', analysis)
        return d

    def _export_sample(self, sample):
        d = {}
        self._fill_id_slots(sample, d)
        d['sample_type'] = sample.data_type
        d['assay_type'] = sample.assay.data_type
        self._add(d, 'strain', sample.strain)
        self._add(d, 'genotype', sample.genotype)
        d['tissue'] = sample.anatomical_part
        d['stage'] = sample.developmental_stage
        self._add(d, 'sex', sample.sex)
        self._add(d, 'technical_reference', sample.assay.technical_reference)
        self._add(d, 'biological_reference', sample.assay.biological_reference)
        if sample.entities:
            entities = []
            for entity, entity_type in sample.entities:
                entities.append({'entity': entity, 'entity_type': entity_type})
            d['experimental_factors'] = entities
        self._add(d, 'collection_fbcv', sample.fbcv)
        self._add(d, 'assay_fbcv', sample.assay.fbcv)
        self._add(d, 'collection_protocol', sample.collection_protocol)
        return d

    def _export_analysis(self, result):
        d = {}
        self._fill_id_slots(result, d)
        d['analysis_type'] = result.data_type
        self._add(d, 'analysis_protocol', result.analysis_protocol)
        d['cell_count'] = result.count
        clusters = []
        for cluster in result.clusters:
            c = {}
            self._fill_id_slots(cluster, c)
            c['cell_type'] = cluster.cell_type
            c['cell_count'] = cluster.count
            clusters.append(c)
        d['clusters'] = clusters
        return d
