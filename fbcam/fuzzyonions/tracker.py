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

    @classmethod
    def values(cls):
        """Gets possible JSON values."""

        return [str(v) for v in cls]

    @classmethod
    def from_click(cls, ctx, param, value):
        """Helper method to get a value from a click Choice."""

        if value:
            return cls.from_str(value)
        else:
            return None


class SceaStatus(JsonEnum):
    """Status of a dataset in the SCEA pipeline."""

    UNKNOWN = 0,
    IN_QUEUE = 1,
    IN_PROGRESS = 2,
    IN_STAGING = 3,
    IN_PRODUCTION = 4


class FlyBaseRecordStatus(JsonEnum):
    """Status of a FlyBase record for a dataset."""

    UNKNOWN = 0,
    IN_PROGRESS = 1,
    READY = 2,
    LOADING = 3,
    LOADED = 4,
    IN_PRODUCTION = 5


class FlyBaseEvaluation(JsonEnum):
    """FlyBase decision about integrating a dataset."""

    UNKNOWN = 0,
    INCLUDE = 1,
    EXCLUDE = 2,
    HOLD = 3


class CorrectionStatus(JsonEnum):
    """Status of a correction set."""

    UNKNOWN = 0,
    NOT_NEEDED = 1,
    TODO = 2,
    DONE = 3


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
        self._flybase = None

    @property
    def scea(self):
        if self._scea is None:
            self._scea = SceaData(basedict=self._data['scea'])
        return self._scea

    @property
    def flybase(self):
        if self._flybase is None:
            self._flybase = FlyBaseData(basedict=self._data.get('flybase', {}))
        return self._flybase

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


class FlyBaseData(object):
    """Represents the FlyBase view of a dataset."""

    def __init__(self, basedict={}):
        self._data = basedict
        self._corrections = None

    @property
    def reference(self):
        return self._data.get('FBrf')

    @property
    def decision(self):
        return FlyBaseEvaluation.from_str(
            self._data.get('evaluation', {}).get('decision', 'unknown'))

    @property
    def is_accepted(self):
        return self.decision == FlyBaseEvaluation.INCLUDE

    @property
    def is_rejected(self):
        return self.decision == FlyBaseEvaluation.EXCLUDE

    @property
    def is_held(self):
        return self.decision == FlyBaseEvaluation.HOLD

    @property
    def comment(self):
        return self._data.get('evaluation', {}).get('comment')

    @property
    def record_name(self):
        return self._data.get('record', {}).get('name')

    @property
    def record_status(self):
        return FlyBaseRecordStatus.from_str(
            self._data.get('record', {}).get('status', 'unknown'))

    @property
    def sumexpr_status(self):
        return FlyBaseRecordStatus.from_str(
            self._data.get('sumexpr', {}).get('status', 'unknown'))

    @property
    def corrections(self):
        if self._corrections is None:
            self._corrections = CorrectionData(
                basedict=self._data.get('corrections', {}))
        return self._corrections


class CorrectionData():
    """Represents the informations about a correction set."""

    def __init__(self, basedict={}):
        self._data = basedict

    @property
    def status(self):
        return CorrectionStatus.from_str(self._data.get('status', 'unknown'))

    @property
    def submitted(self):
        dt = self._data.get('submitted')
        if dt is None:
            return None
        return datetime.strptime(dt, '%Y-%m-%d')

    def to_string(self):
        if self.status == CorrectionStatus.NOT_NEEDED:
            return "no corrections needed"
        elif self.status == CorrectionStatus.TODO:
            return "corrections needed"
        elif self.status == CorrectionStatus.DONE:
            dt = self.submitted
            if dt:
                return f"corrections submitted on {dt:%Y-%m-%d}"
            else:
                return "corrections done, not submitted"
        else:
            return "unknown"


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

        if len(self._data) == 0:
            return "unknown"
        elif self.is_available:
            return "available"
        elif not self.exists:
            return "inexistent"
        else:
            req_date = self.date_requested
            if req_date is None:
                return "to be requested from upstream"
            else:
                return f"requested on {req_date:%Y-%m-%d}"


@click.group(invoke_without_command=True)
@click.pass_context
def tracker(ctx):
    """Track datasets."""

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-tracker> ")
        shell.cmdloop()


@tracker.command('list')
@click.option('--scea-status', type=click.Choice(SceaStatus.values()),
              callback=SceaStatus.from_click,
              help="Filter datasets on their SCEA status.")
@click.option('--fb-decision', type=click.Choice(FlyBaseEvaluation.values()),
              callback=FlyBaseEvaluation.from_click,
              help="Filter datasets on the FlyBase decision status.")
@click.option('--fb-status', type=click.Choice(FlyBaseRecordStatus.values()),
              callback=FlyBaseRecordStatus.from_click,
              help="Filter datasets on their FlyBase curation status.")
@click.pass_obj
def list_tracked_datasets(ctx, scea_status, fb_decision, fb_status):
    """List tracked datasets."""

    for ds in ctx.tracker.datasets:
        if scea_status and ds.scea.status != scea_status:
            continue

        if fb_decision and ds.flybase.decision != fb_decision:
            continue

        if fb_status and ds.flybase.record_status != fb_status:
            continue

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

    print(f"FlyBase evaluation: {ds.flybase.decision}", end='')
    if ds.flybase.comment:
        print(f" ({ds.flybase.comment})", end='')
    print()
    if ds.flybase.decision == FlyBaseEvaluation.INCLUDE:
        print(f"Record: {ds.flybase.record_status}", end='')
        if ds.flybase.record_name:
            print(f" ({ds.flybase.record_name})", end='')
        print()

        if ds.cell_types.is_available:
            print(f"Summarised expression table: {ds.flybase.sumexpr_status}")

    print(f"Corrections: {ds.flybase.corrections.to_string()}")
