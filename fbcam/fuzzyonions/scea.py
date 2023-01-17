# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2021 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

"""Helper module to deal with SCEA datasets."""

from enum import IntFlag
import logging
import os
import json
from datetime import datetime
from zipfile import ZipFile, BadZipFile

import requests
import pandas

from fbcam.fuzzyonions.matrixmarket import MatrixMarketFile


class DataType(IntFlag):
    EXPERIMENT_METADATA = 1
    EXPERIMENT_DESIGN = 2
    CLUSTERS = 4
    MARKER_GENES = 8
    NORMALISED_EXPRESSION_DATA = 16
    RAW_EXPRESSION_DATA = 32


class FileStore(object):
    """Represents a local file-based store for SCEA data."""

    def __init__(self, directory, staging=False):
        self._directory = directory
        self._staging = staging
        self._downloader = Downloader(directory, staging)
        self._datatypes = (
            DataType.EXPERIMENT_DESIGN
            | DataType.EXPERIMENT_METADATA
            | DataType.NORMALISED_EXPRESSION_DATA
            | DataType.RAW_EXPRESSION_DATA
        )

        if not os.path.exists(directory):
            os.makedirs(directory)

        self._datasets = None

    @property
    def datasets(self):
        """Gets all the datasets in the store."""

        if self._datasets is None:
            self._refresh()
        return [d for d in self._datasets.values()]

    def get(self, dsid, download=True):
        """Gets a single dataset.

        :param dsid: the dataset ID
        :param download: if True (the default), the dataset will be
            downloaded if it is not already available locally
        """

        if self._datasets is None:
            self._refresh()

        if dsid not in self._datasets.keys():
            if not download:
                return None
            elif not self._downloader.get(dsid, self._datatypes):
                return None
            self._datasets[dsid] = Dataset(
                os.path.join(self._directory, dsid), self._staging
            )

        return self._datasets[dsid]

    def delete(self, dsids):
        """Remove the specified datasets from the local file store."""

        if self._datasets is None:
            self._refresh()

        for dsid in dsids:
            if dsid in self._datasets:
                self._datasets[dsid].delete()
        self._refresh()

    def _refresh(self):
        self._datasets = {}
        for dsid in os.listdir(self._directory):
            fullpath = os.path.join(self._directory, dsid)
            if os.path.isdir(fullpath):
                self._datasets[dsid] = Dataset(
                    os.path.join(self._directory, dsid), self._staging
                )

    def get_experiments_list(self, force_refresh=False):
        expfile = os.path.join(self._directory, 'experiments.json')
        mtime = None

        if force_refresh or not os.path.exists(expfile):
            refresh = True
        else:
            mtime = os.path.getmtime(expfile)
            refresh = mtime + 86400 < datetime.now().timestamp()

        if refresh:
            if mtime is not None:
                mtime = datetime.fromtimestamp(mtime).astimezone()
            self._downloader.get_experiments_list(since=mtime)

        with open(expfile, 'r') as fd:
            return json.load(fd)


class Dataset(object):
    """Represents a single SCEA dataset."""

    FILES = {
        DataType.EXPERIMENT_METADATA: '{}.idf.txt',
        DataType.EXPERIMENT_DESIGN: 'experiment-design.tsv',
        DataType.CLUSTERS: '{}.clusters.tsv',
        DataType.MARKER_GENES: '{}-marker-genes-files',
        DataType.NORMALISED_EXPRESSION_DATA: '{}.aggregated_filtered_normalised_counts.mtx',
        DataType.RAW_EXPRESSION_DATA: '{}.aggregated_filtered_counts.mtx',
    }

    def __init__(self, directory, staging=False):
        self._directory = directory
        self._id = os.path.basename(directory)
        self._files = os.listdir(directory)
        self._expdesign = None
        self._staging = staging

    @property
    def id(self):
        return self._id

    @property
    def staging(self):
        return self._staging

    @property
    def available_data(self):
        """Gets a list of locally-available data for this dataset."""

        return [
            d for d, f in Dataset.FILES.items() if f.format(self._id) in self._files
        ]

    @property
    def experiment_design(self):
        """Gets the "experiment design" table of the dataset."""

        if self._expdesign is None:
            f = Dataset.FILES[DataType.EXPERIMENT_DESIGN].format(self._id)
            self._expdesign = pandas.read_csv(
                os.path.join(self._directory, f), sep='\t', low_memory=False
            )
        return self._expdesign

    def get_expression_matrix(self, raw=True):
        """Gets the expression matrix table, as a file."""

        return MatrixMarketFile(self._get_expression_matrix_fullname(raw))

    def delete(self):
        """Delete this dataset from its file store."""

        for file in self._files:
            os.remove(os.path.join(self._directory, file))
        os.rmdir(self._directory)

    def _get_expression_matrix_fullname(self, raw=True):
        """Gets the full pathname to the expression matrix file."""

        dt = DataType.RAW_EXPRESSION_DATA
        if not raw:
            dt = DataType.NORMALISED_EXPRESSION_DATA

        basename = Dataset.FILES[dt].format(self._id)
        fullname = os.path.join(self._directory, basename)

        return fullname


class Downloader(object):
    """Helper object to download data from the SCEA."""

    BASE_URL = 'https://www.ebi.ac.uk/gxa/sc/'
    BASE_URL_DEV = 'https://wwwdev.ebi.ac.uk/gxa/sc/'
    FILES = {
        DataType.EXPERIMENT_METADATA: ('experiment-metadata', True),
        DataType.EXPERIMENT_DESIGN: ('experiment-design', False),
        DataType.CLUSTERS: ('cluster', False),
        DataType.MARKER_GENES: ('marker-genes', True),
        DataType.NORMALISED_EXPRESSION_DATA: ('normalised', True),
        DataType.RAW_EXPRESSION_DATA: ('quantification-raw', True),
    }

    def __init__(self, destination, staging=False):
        """Creates a new Downloader object.

        :param destination: the directory to store the downloaded
            files into
        :param base: an alternative URL to download from
        """

        self._outdir = destination
        if staging:
            self._baseurl = self.BASE_URL_DEV
        else:
            self._baseurl = self.BASE_URL

    def get(self, dsid, data_type):
        """Download a dataset file from the SCEA.

        :param dsid: the dataset ID
        :param data_type: the type of data to download
        :return: True if the download was successful, False otherwise
        """

        errors = 0

        if not self._check_dataset(dsid):
            return False

        for value in DataType:
            if data_type & value:
                try:
                    filename, is_zip = Downloader.FILES[value]
                    self._download(dsid, filename, is_zip)
                except requests.RequestException as e:
                    logging.info(f"Cannot download {filename}: {e}")
                    errors += 1
                except BadZipFile as e:
                    logging.info(f"Cannot unzip {filename}: {e}")
                    errors += 1

        return errors == 0

    def get_experiments_list(self, since=None):
        """Download the list of experiments, as a JSON file.

        :param since: a DateTime object; only download the file if the
            version on the server is more recent than the indicate date
        :return: True if the download was successful (or there was no
            newer version on the server), False otherwise
        """

        url = f'{self._baseurl}json/experiments'
        dest = f'{self._outdir}/experiments.json'
        headers = {}
        if since is not None:
            headers['If-Modified-Since'] = since.strftime('%a, %d %b %Y %H:%M:%S %Z')
        try:
            with requests.get(url, headers=headers) as response:
                if response.status_code == 304:
                    return True
                data = [
                    e
                    for e in json.loads(response.content)['experiments']
                    if e['species'] == 'Drosophila melanogaster'
                ]
                with open(dest, 'w') as fd:
                    json.dump(data, fd, indent=2)
        except requests.RequestException as e:
            logging.info(f"Cannot download experiments list: {e}")
            return False

        return True

    def _check_dataset(self, dsid):
        url = f'{self._baseurl}experiments/{dsid}/results'
        try:
            with requests.head(url) as response:
                return response.ok
        except requests.RequestException:
            return False

    def _download(self, dsid, file_type, is_zip=False):
        z = 'zip' if is_zip else ''
        ext = 'zip' if is_zip else 'tsv'
        url = f'{self._baseurl}experiment/{dsid}/download/{z}?fileType={file_type}'
        destdir = f'{self._outdir}/{dsid}'
        dest = f'{destdir}/{file_type}.{ext}'

        if not os.path.exists(destdir):
            os.mkdir(destdir)

        logging.info(f"Downloading {url}")

        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(dest, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=None):
                fd.write(chunk)

        if is_zip:
            archive = ZipFile(dest)
            for member in archive.namelist():
                archive.extract(member, path=destdir)
            os.unlink(dest, dir_fd=None)
