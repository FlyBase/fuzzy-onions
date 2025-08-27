# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2025 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

from pronto import Ontology

PURL_BASE = 'http://purl.obolibrary.org/obo/'


class OntologyStore(object):
    """Centralised access to all ontologies."""

    def __init__(self, config):
        self.fbbt = self._get_ontology('fbbt', config)
        self.fbdv = self._get_ontology('fbdv', config)
        self.fbcv = self._get_ontology('fbcv', config)

    def _get_ontology(self, name, config):
        path = config.get('ontologies', name, fallback=f'{PURL_BASE}{name}.obo')
        return CachedOntology(path)


class CachedOntology(object):
    """Helper object to access an ontology."""

    def __init__(self, path):
        """Creates a new instance.

        :param path: The path to the file containing the ontology, or a URL
                     pointing to the online location of the ontology.
        """
        self._path = path
        self._backend = None
        self._ids_by_label = None
        self._labels_by_exact_syn = None

    @property
    def backend(self):
        """The Pronto ontology object."""
        if self._backend is None:
            self._backend = Ontology(self._path)
        return self._backend

    @property
    def ids(self):
        """A dictionary associating labels to IDs."""
        if self._ids_by_label is None:
            self._ids_by_label = {}
            for term in self.backend.terms():
                self._ids_by_label[term.name] = term.id
        return self._ids_by_label

    @property
    def synonyms(self):
        """A dictionary associating exact synonyms to canonical labels."""
        if self._labels_by_exact_syn is None:
            self._labels_by_exact_syn = {}
            for term in self.backend.terms():
                for synonym in [s for s in term.synonyms if s.scope == 'EXACT']:
                    self._labels_by_exact_syn[synonym.description] = term.name
        return self._labels_by_exact_syn
