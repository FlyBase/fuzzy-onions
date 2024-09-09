# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2021 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import os
from enum import Enum


class ProformaType(Enum):
    PUB_MINI = (0,)
    DATASET = 1


class ProformaGeneratorBuilder(object):

    _PROFORMAE = {
        ProformaType.PUB_MINI: ('pub_mini', 'PUBLICATION'),
        ProformaType.DATASET: ('dataset_master', 'DATASET/COLLECTION'),
    }

    def __init__(self, directory, output):
        self._directory = directory
        self._output = output
        self._generators = {}

    def get_generator(self, proforma_type=ProformaType.DATASET):

        if not proforma_type in self._generators:
            filename, name = self._PROFORMAE[proforma_type]
            pathname = os.path.join(self._directory, filename + '.pro')

            generator = ProformaGenerator(pathname, name, self._output)
            self._generators[proforma_type] = generator

        generator = self._generators[proforma_type]

        return generator


class ProformaGenerator(object):
    def __init__(self, source, name, output):
        self._source = source
        self._name = name
        self._fields = None
        self._fp = output

    @property
    def fields(self):
        if not self._fields:
            self._fields = self._parse_proforma(self._source, self._name)
        return self._fields

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
        self._fp.write(
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
        )

    def write_terminator(self):
        self._fp.write(
            '!!!!!!!!!!!!!!!!!! END OF RECORD FOR THIS PUBLICATION !!!!!!!!!!!!!!!!!!!!\n'
        )

    def write_field(self, fid, value=None):
        if self._is_empty(value):
            # Do not output the field if it is explicitly empty
            return

        self._fp.write(self.fields[fid])
        if value is not None:
            if isinstance(value, list):
                self._fp.write('\n'.join(value))
            else:
                self._fp.write(str(value))
        self._fp.write('\n')

    def _is_empty(self, value):
        if value is None:
            return False
        elif isinstance(value, list):
            return len(value) == 0
        elif isinstance(value, str):
            return len(value) == 0
        else:
            return False
