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


class CellTypeAvailability(JsonEnum):
    """Status of cell type annotations."""

    UNKNOWN = 0,
    INEXISTENT = 1,
    INPUT_ONLY = 2,
    UPSTREAM = 3,
    OBTAINED = 4,
    AVAILABLE = 5


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
                    self._datasets.append(TrackedDataset.from_dict(d))

    def get_dataset(self, dsid):
        """Gets a dataset from its SCEA ID."""

        for ds in self.datasets:
            if ds.scea.identifier == dsid:
                return ds
        return None

    def save(self, newfile=None):
        """Writes the tracking database to a file."""

        if newfile is None:
            newfile = self._db_file
        with open(newfile, 'w') as f:
            json.dump([d.to_dict() for d in self.datasets], f, indent=2)

    def add_dataset(self, dsid, staging=False):
        """Adds a dataset to track."""

        if self.get_dataset(dsid):
            return

        status = SceaStatus.IN_PRODUCTION
        if staging:
            status = SceaStatus.IN_STAGING

        ds = { 'scea': { 'dataset_id': dsid, 'status': str(status)} }
        self._datasets.append(TrackedDataset.from_dict(ds))


class TrackedDataset(object):
    """Represents a single dataset."""

    def __init__(self):
        self._scea = SceaData()
        self._flybase = None
        self._ctypes = None

    @property
    def scea(self):
        return self._scea

    @property
    def flybase(self):
        if self._flybase is None:
            self._flybase = FlyBaseData()
        return self._flybase

    @property
    def cell_types(self):
        if self._ctypes is None:
            self._ctypes = CellTypesData()
        return self._ctypes

    def to_dict(self):
        d = { 'scea': self._scea.to_dict() }
        if not self.flybase.is_empty:
            d['flybase'] = self.flybase.to_dict()
        if self.cell_types._status != CellTypeAvailability.UNKNOWN:
            d['cell_types_data'] = self.cell_types.to_dict()
        return d

    @classmethod
    def from_dict(cls, data):
        new = cls()

        new._scea = SceaData.from_dict(data['scea'])
        if 'flybase' in data:
            new._flybase = FlyBaseData.from_dict(data['flybase'])
        if 'cell_types_data' in data:
            new._ctypes = CellTypesData.from_dict(data['cell_types_data'])

        return new


class SceaData(object):
    """Represents the SCEA view of a dataset."""

    def __init__(self):
        self._dataset_id = '<unknown id>'
        self._upstream_id = None
        self._status = SceaStatus.UNKNOWN

    @property
    def identifier(self):
        return self._dataset_id

    @property
    def status(self):
        return self._status

    @property
    def upstream_identifier(self):
        return self._upstream_id

    def to_dict(self):
        d = { 'dataset_id': self._dataset_id }
        if self._upstream_id:
            d['upstream_id'] = self._upstream_id
        if self._status != SceaStatus.UNKNOWN:
            d['status'] = str(self._status)
        return d

    @classmethod
    def from_dict(cls, data):
        new = cls()
        new._dataset_id = data.get('dataset_id', '<unknown id>')
        new._upstream_id = data.get('upstream_id')
        new._status = SceaStatus.from_str(data.get('status', 'unknown'))
        return new


class FlyBaseData(object):
    """Represents the FlyBase view of a dataset."""

    def __init__(self):
        self._reference = None
        self._decision = FlyBaseEvaluation.UNKNOWN
        self._comment = None
        self._name = None
        self._status = FlyBaseRecordStatus.UNKNOWN
        self._sumexpr = FlyBaseRecordStatus.UNKNOWN
        self._corrections = CorrectionData()

    @property
    def is_empty(self):
        return (self._reference is None and
                self._decision == FlyBaseEvaluation.UNKNOWN and
                self._status == FlyBaseRecordStatus.UNKNOWN and
                self._sumexpr == FlyBaseRecordStatus.UNKNOWN)

    @property
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, value):
        self._reference = value

    @property
    def decision(self):
        return self._decision

    def decide(self, decision, comment=None):
        self._decision = decision
        self._comment = comment

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
        return self._comment

    @property
    def record_name(self):
        return self._name

    @record_name.setter
    def record_name(self, value):
        self._name = value

    @property
    def record_status(self):
        return self._status

    def update_record(self, status, name=None):
        self._status = status
        if name:
            self._name = name

    @property
    def sumexpr_status(self):
        return self._sumexpr

    @property
    def corrections(self):
        return self._corrections

    def get_status_string(self):
        if self.decision == FlyBaseEvaluation.UNKNOWN:
            return "undecided"
        elif self.decision in [FlyBaseEvaluation.EXCLUDE,
                                FlyBaseEvaluation.HOLD]:
            return str(self.decision)
        elif self.record_status == FlyBaseRecordStatus.UNKNOWN:
            return "todo"
        else:
            return str(self.record_status)

    def to_dict(self):
        d = {}
        if self._reference:
            d['fbrf'] = self._reference
        if self._decision != FlyBaseEvaluation.UNKNOWN:
            d['evaluation'] = { 'decision': str(self._decision) }
            if self._comment:
                d['evaluation']['comment'] = self._comment
        if self._status != FlyBaseRecordStatus.UNKNOWN:
            d['record'] = { 'status': str(self._status) }
            if self._name:
                d['record']['name'] = self._name
        if self._sumexpr != FlyBaseRecordStatus.UNKNOWN:
            d['sumexpr'] = { 'status': str(self._sumexpr)}
        if self.corrections.status != CorrectionStatus.UNKNOWN:
            d['corrections'] = self.corrections.to_dict()
        return d

    @classmethod
    def from_dict(cls, data):
        new = cls()

        new._reference = data.get('fbrf')

        if 'evaluation' in data:
            new._decision = FlyBaseEvaluation.from_str(
                data['evaluation'].get('decision', 'unknown'))
            new._comment = data['evaluation'].get('comment')

        if 'record' in data:
            new._status = FlyBaseRecordStatus.from_str(
                data['record'].get('status', 'unknown'))
            new._name = data['record'].get('name')

        if 'sumexpr' in data:
            new._sumexpr = FlyBaseRecordStatus.from_str(
                data['sumexpr'].get('status', 'unknown'))

        if 'corrections' in data:
            new._corrections = CorrectionData.from_dict(data['corrections'])

        return new


class CorrectionData():
    """Represents the informations about a correction set."""

    def __init__(self):
        self._status = CorrectionStatus.UNKNOWN
        self._submitted = None

    @property
    def status(self):
        return self._status

    @property
    def submitted(self):
        return self._submitted

    def to_string(self):
        if self.status == CorrectionStatus.NOT_NEEDED:
            return "none needed"
        elif self.status == CorrectionStatus.TODO:
            return "todo"
        elif self.status == CorrectionStatus.DONE:
            dt = self.submitted
            if dt:
                return f"submitted on {dt:%Y-%m-%d}"
            else:
                return "unsubmitted"
        else:
            return "unknown"

    def to_dict(self):
        d = { 'status': str(self._status) }
        if self._submitted:
            d['submitted'] = self._submitted.strftime('%Y-%m-%d')
        return d

    @classmethod
    def from_dict(cls, data):
        new = cls()
        new._status = CorrectionStatus.from_str(data.get('status', 'unknown'))
        if 'submitted' in data:
            new._submitted = datetime.strptime(data['submitted'], '%Y-%m-%d')
        return new


class CellTypesData(object):
    """Represents the informations on cell types annotations."""

    def __init__(self):
        self._status = CellTypeAvailability.UNKNOWN
        self._requested = None

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, value):
        self._status = value

    @property
    def date_requested(self):
        """When have annotations been requested from upstream?"""

        return self._requested

    @property
    def need_request(self):
        return (self._status == CellTypeAvailability.UPSTREAM and
                self._requested is None)

    def set_to_request(self):
        """Marks that cell types are to be requested."""

        self._status = CellTypeAvailability.UPSTREAM
        self._requested = None

    def set_requested(self, date=datetime.today()):
        """Marks that cell types have been requested."""

        self._status = CellTypeAvailability.UPSTREAM
        self._requested = date

    def to_string(self):
        """Gets a one-line human-readable representation."""

        if self._status == CellTypeAvailability.UNKNOWN:
            return "unknown"
        elif self._status == CellTypeAvailability.AVAILABLE:
            return "available"
        elif self._status == CellTypeAvailability.INEXISTENT:
            return "inexistent"
        elif self._status == CellTypeAvailability.INPUT_ONLY:
            return "input types only"
        elif self._requested is None:
            return "to be requested from upstream"
        else:
            return f"requested on {self._requested:%Y-%m-%d}"

    def to_dict(self):
        d = { 'status': str(self._status) }
        if self._requested:
            d['requested'] = self._requested.strftime('%Y-%m-%d')
        return d

    @classmethod
    def from_dict(cls, data):
        new = cls()
        new._status = CellTypeAvailability.from_str(data.get('status', 'unknown'))
        if 'requested' in data:
            new._requested = datetime.strptime(data['requested'], '%Y-%m-%d')

        return new


@click.group(invoke_without_command=True)
@click.pass_context
def tracker(ctx):
    """Track datasets."""

    if not ctx.invoked_subcommand:
        ctx.obj.in_tracker_shell = True
        shell = make_click_shell(ctx, prompt="fzo-tracker> ")
        shell.cmdloop()
    else:
        ctx.obj.in_tracker_shell = False


@tracker.command('list')
@click.option('--long', '-l', is_flag=True, default=False,
              help="Show more details than just the dataset ID.")
@click.option('--tab', '-t', is_flag=True, default=False,
              help="Use tab-separated columns.")
@click.option('--scea-status', type=click.Choice(SceaStatus.values()),
              callback=SceaStatus.from_click,
              help="Filter datasets on their SCEA status.")
@click.option('--fb-decision', type=click.Choice(FlyBaseEvaluation.values()),
              callback=FlyBaseEvaluation.from_click,
              help="Filter datasets on the FlyBase decision status.")
@click.option('--fb-status', type=click.Choice(FlyBaseRecordStatus.values()),
              callback=FlyBaseRecordStatus.from_click,
              help="Filter datasets on their FlyBase curation status.")
@click.option('--ct-status', type=click.Choice(CellTypeAvailability.values()),
              callback=CellTypeAvailability.from_click,
              help="Filter datasets on cell types availability.")
@click.pass_obj
def list_tracked_datasets(ctx, long, tab, scea_status, fb_decision, fb_status,
                          ct_status):
    """List tracked datasets."""

    if long:
        print("Dataset ID       EBI Status       FlyBase Status   Cell Types                       Corrections")
    elif tab:
        print("Dataset ID\tEBI Status\tFlyBase Status\tCell Types\tCorrections")

    for ds in ctx.tracker.datasets:
        if scea_status and ds.scea.status != scea_status:
            continue

        if fb_decision and ds.flybase.decision != fb_decision:
            continue

        if fb_status and ds.flybase.record_status != fb_status:
            continue

        if ct_status and ds.cell_types.status != ct_status:
            continue

        if long:
            print(f"{ds.scea.identifier:16} {ds.scea.status:16} {ds.flybase.get_status_string():16} {ds.cell_types.to_string():32} {ds.flybase.corrections.to_string()}")
        elif tab:
            print(f"{ds.scea.identifier}\t{ds.scea.status}\t{ds.flybase.get_status_string()}\t{ds.cell_types.to_string()}\t{ds.flybase.corrections.to_string()}")
        else:
            print(ds.scea.identifier)


@tracker.command()
@click.argument('dsid')
@click.pass_obj
def show(ctx, dsid):
    """Show information about a tracked dataset."""

    ds = ctx.tracker.get_dataset(dsid)
    if not ds:
        raise click.ClickException("Invalid Dataset ID.")

    print(f"Dataset ID: {ds.scea.identifier}")
    print(f"SCEA status: {ds.scea.status}")
    print(f"Cell type annotations: {ds.cell_types.to_string()}")

    print(f"FlyBase reference: {ds.flybase.reference}")
    print(f"FlyBase evaluation: {ds.flybase.decision}", end='')
    if ds.flybase.comment:
        print(f" ({ds.flybase.comment})", end='')
    print()
    if ds.flybase.decision == FlyBaseEvaluation.INCLUDE:
        print(f"Record: {ds.flybase.record_status}", end='')
        if ds.flybase.record_name:
            print(f" ({ds.flybase.record_name})", end='')
        print()

        if ds.cell_types.status == CellTypeAvailability.AVAILABLE:
            print(f"Summarised expression table: {ds.flybase.sumexpr_status}")

    print(f"Corrections: {ds.flybase.corrections.to_string()}")


@tracker.command()
@click.argument('dsid')
@click.pass_obj
def accept(ctx, dsid):
    """Mark a dataset as accepted for curation in FlyBase."""

    ds = ctx.tracker.get_dataset(dsid)
    if not ds:
        raise click.ClickException("Invalid Dataset ID.")

    ds.flybase.decide(FlyBaseEvaluation.INCLUDE)
    if not ctx.in_tracker_shell:
        ctx.tracker.save()


@tracker.command()
@click.pass_obj
def add(ctx, dsid):
    """Add a Dataset ID to track."""

    ctx.tracker.add_dataset(dsid)


@tracker.command()
@click.argument('dsid')
@click.option('--cell-types',
              type=click.Choice(CellTypeAvailability.values()),
              callback=CellTypeAvailability.from_click,
              help="Update what’s known about cell type annotations.")
@click.option('--decide', 'decision',
              type=click.Choice(FlyBaseEvaluation.values()),
              callback=FlyBaseEvaluation.from_click,
              help="Set the FlyBase decision status.")
@click.option('--comment', help="Set the comment from the FlyBase curator.")
@click.option('--reference', help="Set the associated FlyBase reference.")
@click.option('--record', help="Set the name of the FlyBase record.")
@click.pass_obj
def update(ctx, dsid, cell_types, decision, comment, reference, record):
    """Update informations about a dataset."""

    ds = ctx.tracker.get_dataset(dsid)
    if not ds:
        raise click.ClickException("Invalid Dataset ID.")

    if cell_types:
        ds.cell_types.status = cell_types

    if comment is not None and decision is None:
        decision = FlyBaseEvaluation.UNKNOWN
    if decision:
        ds.flybase.decide(decision, comment)

    if reference:
        ds.flybase.reference = reference
    if record:
        ds.flybase.record_name = record

    if not ctx.in_tracker_shell:
        ctx.tracker.save()


@tracker.command()
@click.option('--filename', '-f',
              type=click.Path(exists=False, writable=True, dir_okay=False),
              help="Write to a different file.")
@click.pass_obj
def save(ctx, filename):
    """Write modification to the tracking data."""

    ctx.tracker.save(filename)
