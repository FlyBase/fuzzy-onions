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

import os
from enum import Enum


class ProformaType(Enum):
    PUB_MINI = 0,
    DATASET = 1


TEMPLATES = {
    'pub_mini': {
        'proforma': ProformaType.PUB_MINI,
        'contents': [
            ('P22', 'TODO: FBrfXXXXXXX'),
            ('P2', 'TODO: <Journal abbreviation>'),
            ('P19', None)
            ]
        },
    'dataset/project': {
        'proforma': ProformaType.DATASET,
        'contents': [
            ('LC1f', 'new'),
            ('LC1a', None),
            ('LC6g', 'TODO: Single-cell RNA-seq study of <tissue and stage> [upon condition]'),
            ('LC6a', None),
            ('LC2a', 'project ; FBcv:0003023'),
            ('LC2b', 'transcriptome ; FBcv:0003034'),
            ('LC13a', None),
            ('LC13b', None),
            ('LC13c', 'TODO: [GO terms for biological processes studied]'),
            ('LC13d', None),
            ('LC4a', 'Dmel'),
            ('LC4i', None),
            ('LC6d', 'N'),
            ('LC11m', 'TODO: <FBcv terms from \'study design\' (FBcv:0003130)>'),
            ('LC6b', 'TODO: <How the biological material was processed>'),
            ('LC11c', 'TODO: <How the mRNA libraries were sequenced>'),
            ('LC11e', 'TODO: <How the sequencing data were analysed>'),
            ('LC7a', 'The EMBL-EBI\'s Single Cell Expression Atlas provides cell-level annotations, clustering data, raw and normalised read counts, and putative marker genes.'),
            ('LC99a', None),
            ('LC99b', 'EMBL-EBI Single Cell Expression Atlas Datasets'),
            ('LC8c', 'TODO: [Name of the research group](URL of their website)'),
            ('LC10', None)
            ]
        },
    'dataset/biosample': {
        'proforma': ProformaType.DATASET,
        'contents': [
            ('LC1f', 'new'),
            ('LC1a', None),
            ('LC6g', None),
            ('LC2a', 'biosample ; FBcv:0003024'),
            ('LC2b', None),
            ('LC3', None),
            ('LC4a', 'Dmel'),
            ('LC4h', 'TODO: [as needed]'),
            ('LC4f', 'TODO: [as needed]'),
            ('LC4g', None),
            ('LC4i', None),
            ('LC12a', 'TODO: [as needed]'),
            ('LC12b', 'TODO: [as needed]'),
            ('LC6d', 'N'),
            ('LC11m', None),
            ('LC11a', None),
            ('LC10', None)
            ]
        },
    'dataset/assay': {
        'proforma': ProformaType.DATASET,
        'contents': [
            ('LC1f', 'new'),
            ('LC1a', None),
            ('LC6g', None),
            ('LC2a', 'assay ; FBcv:0003025'),
            ('LC2b', None),
            ('LC3', None),
            ('LC14a', None),
            ('LC14d', 'TODO: [LC1a symbol of technical reference assay, if relevant]'),
            ('LC14e', 'TODO: [LC1a symbol of biological reference assay, if relevant]'),
            ('LC4a', 'Dmel'),
            ('LC6d', 'N'),
            ('LC6e', None),
            ('LC6f', None),
            ('LC11m', 'TODO: <FBcv terms from \'assay method\' (FBcv:0003208)>'),
            ('LC10', None)
            ]
        },
    'dataset/result': {
        'proforma': ProformaType.DATASET,
        'contents': [
            ('LC1f', 'new'),
            ('LC1a', None),
            ('LC6g', None),
            ('LC2a', 'result ; FBcv:0003026'),
            ('LC2b', None),
            ('LC3', None),
            ('LC14b', None),
            ('LC14h', 'Dmel R6.32'),
            ('LC4a', 'Dmel'),
            ('LC6d', 'N'),
            ('LC6e', None),
            ('LC6f', None),
            ('LC10', None)
            ]
        },
    'dataset/subresult': {
        'proforma': ProformaType.DATASET,
        'contents': [
            ('LC1f', 'new'),
            ('LC1a', None),
            ('LC6g', None),
            ('LC2a', 'result ; FBcv:0003026'),
            ('LC2b', None),
            ('LC3', None),
            ('LC4a', 'Dmel'),
            ('LC4g', None),
            ('LC6d', 'N'),
            ('LC6e', None),
            ('LC6f', None),
            ('LC10', None)
            ]
        }
    }


class ProformaGeneratorBuilder(object):

    _PROFORMAE = {
        ProformaType.PUB_MINI: ('pub_mini', 'PUBLICATION'),
        ProformaType.DATASET: ('dataset_master', 'DATASET/COLLECTION')
        }

    def __init__(self, directory, output):
        self._directory = directory
        self._output = output
        self._generators = {}

    def get_generator(self, proforma_type=None, template=None):
        if template is not None:
            proforma_type = TEMPLATES[template]['proforma']
            template = TEMPLATES[template]['contents']

        if not proforma_type in self._generators:
            filename, name = self._PROFORMAE[proforma_type]
            pathname = os.path.join(self._directory, filename + '.pro')

            generator = ProformaGenerator(pathname, name, self._output)
            self._generators[proforma_type] = generator

        generator = self._generators[proforma_type]
        if template is not None:
            generator.set_template(template)

        return generator


class ProformaGenerator(object):

    def __init__(self, source, name, output, template=None):
        self._source = source
        self._name = name
        self._fields = None
        self._fp = output
        self._template = template

    @property
    def fields(self):
        if not self._fields:
            self._fields = self._parse_proforma(self._source, self._name)
        return self._fields

    def set_template(self, template):
        self._template = template

    def _parse_proforma(self, filename, name):
        res = {}
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(f'! {name} PROFORMA'):
                    res['header'] = line
                elif line.startswith('! '):
                    fid = line.split('.')[0][2:]
                    res[fid] = line
        return res

    def write_header(self):
        self.write_field('header')
        self.write_blank()

    def write_blank(self):
        self._fp.write('!\n')

    def write_separator(self):
        self._fp.write('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')

    def write_terminator(self):
        self._fp.write('!!!!!!!!!!!!!!!!!! END OF RECORD FOR THIS PUBLICATION !!!!!!!!!!!!!!!!!!!!\n')

    def write_field(self, fid, value=None):
        self._fp.write(self.fields[fid])
        if value is not None:
            self._fp.write(str(value))
        self._fp.write('\n')

    def fill_template(self, values={}, template=None):
        if template is None:
            template = self._template

        self.write_header()

        for field, value in template:
            if field in values:
                if values[field] != '<DEFAULT>':
                    value = values[field]
            self.write_field(field, value)

        self.write_separator()
