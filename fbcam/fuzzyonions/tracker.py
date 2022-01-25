# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
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

from datetime import datetime
from enum import Enum
from os.path import exists
import json

import click
from click_shell import make_click_shell


class JsonEnum(Enum):
    """Helper class to deal with enums in JSON data."""

    def __str__(self):
        """Gets a JSON representation of an enum value."""

        return self.name.lower().replace('_', '-')

    @classmethod
    def from_str(cls, name):
        """Gets a enum value from its JSON representation."""

        return cls[name.replace('-', '_').upper()]


class SceaStatus(JsonEnum):
    """Status of a dataset in the SCEA pipeline."""

    IN_QUEUE = 0,
    IN_PROGRESS = 1,
    IN_STAGING = 2,
    IN_PRODUCTION = 3


class FlyBaseRecordStatus(JsonEnum):
    """Status of a FlyBase record for a dataset."""

    IN_PROGRESS = 0,
    READY = 1,
    LOADING = 2,
    LOADED = 3,
    IN_PRODUCTION = 4


class DatasetTracker(object):
    """Helper object to track datasets.
    
    This class encapsulates a JSON file containing the informations
    to be tracked about the datasets.
    """

    def __init__(self, database):
        self._db_file = database
        self._datasets = None

    @property
    def datasets(self):
        """Gets all the tracked datasets."""

        if self._datasets is None:
            self._load()
        return self._datasets

    def _load(self):
        if not exists(self._db_file):
            self._datasets = []
        else:
            with open(self._db_file, 'r') as f:
                self._datasets = []
                for d in json.load(f):
                    self._datasets.append(TrackedDataset(basedict=d))

    def get_dataset(self, dsid):
        """Gets a dataset from its SCEA ID."""

        for ds in self.datasets:
            if ds.scea.identifier == dsid:
                return ds
        return None


class TrackedDataset(object):
    """Represents a single dataset."""

    def __init__(self, basedict={}):
        self._data = basedict
        self._scea = None
        self._ctypes = None

    @property
    def scea(self):
        if self._scea is None:
            self._scea = SceaData(basedict=self._data['scea'])
        return self._scea

    @property
    def cell_types(self):
        if self._ctypes is None:
            self._ctypes = CellTypesData(basedict=self._data.get('cell_types_data', {}))
        return self._ctypes


class SceaData(object):
    """Represents the SCEA view of a dataset."""

    def __init__(self, basedict={}):
        self._data = basedict

    @property
    def identifier(self):
        return self._data['dataset_id']

    @property
    def status(self):
        return SceaStatus.from_str(self._data['status'])

    @property
    def upstream_identifier(self):
        return self._data.get('upstream_id')


class CellTypesData(object):
    """Represents the informations on cell types annotations."""

    def __init__(self, basedict={}):
        self._data = basedict

    @property
    def exists(self):
        """Have cell types been inferred by the authors?"""

        return self._data.get('exist', 'no') == 'yes'

    @property
    def is_available(self):
        """Are inferred cell types available at the SCEA?"""

        return self._data.get('available', 'no') == 'yes'

    @property
    def date_requested(self):
        """When have annotations been requested from upstream?"""

        dt = self._data.get('requested')
        if dt is None:
            return None
        return datetime.strptime(dt, '%Y-%m-%d')

    def to_string(self):
        """Gets a one-line human-readable representation."""

        if self.is_available:
            return "Available"
        elif not self.exists:
            return "Inexistent"
        else:
            req_date = self.date_requested
            if req_date is None:
                return "To be requested from upstream"
            else:
                return f"Requested on {req_date:%Y-%m-%d}"


@click.group(invoke_without_command=True)
@click.pass_context
def tracker(ctx):
    """Track datasets."""

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-tracker> ")
        shell.cmdloop()


@tracker.command('list')
@click.pass_obj
def list_tracked_datasets(ctx):
    """List tracked datasets."""

    for ds in ctx.tracker.datasets:
        print(ds.scea.identifier)


@tracker.command()
@click.argument('dsid')
@click.pass_obj
def show(ctx, dsid):
    """Show information about a tracked dataset."""

    ds = ctx.tracker.get_dataset(dsid)
    if not ds:
        print(f"Dataset {dsid} not found.")

    print(f"Dataset ID: {ds.scea.identifier}")
    print(f"SCEA status: {ds.scea.status}")

    print(f"Cell type annotations: {ds.cell_types.to_string()}")
