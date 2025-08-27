# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2025 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

from pronto import Ontology


class CachedOntology(object):
    """Helper object to access an ontology.

    FIXME: This is currently using Pronto, maybe consider using OAK instead.
    """

    def __init__(self, path):
        self._path = path
        self._backend = None
        self._reverse_dict = None
        self._synonyms_dict = None

    @property
    def backend(self):
        """The Pronto ontology object."""
        if self._backend is None:
            self._backend = Ontology(self._path)
        return self._backend

    @property
    def ids(self):
        """A dictionary associated labels to IDs."""
        if self._reverse_dict is None:
            self._reverse_dict = {}
            for term in self.backend.terms():
                self._reverse_dict[term.name] = term.id
        return self._reverse_dict
