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


class MatrixMarketFile(object):
    """Helper object to read a MatrixMarket file."""

    def __init__(self, filename):
        self._filename = filename
        self._handle = None
        self._columns = None
        self._rows = None
        self._n_cells = 0
        self._progress_cb = None

    def _get_columns(self):
        with open(f'{self._filename}_cols', 'r') as f:
            return [line.rstrip() for line in f]

    def _get_rows(self):
        with open(f'{self._filename}_rows', 'r') as f:
            return [line.split()[0] for line in f]

    @property
    def rows(self):
        if not self._rows:
            self._rows = self._get_rows()
        return self._rows

    @property
    def columns(self):
        if not self._columns:
            self._columns = self._get_columns()
        return self._columns

    def set_progress_callback(self, callback):
        """Sets a callback function to monitor progress.
        
        Sets a function that will be called for each percent of
        the table read. The function should take a single value
        that will be the percentage of progress.
        """

        self._progress_cb = callback

    def open(self):
        """Opens the MatrixMarket file and check its header."""

        if self._handle is not None:
            return

        self._handle = open(self._filename, 'r')
        header = self._handle.readline()
        if not header.startswith('%%MatrixMarket'):
            raise Exception("Invalid MatrixMarket header.")

        comment = True
        while comment:
            line = self._handle.readline()
            if line[0] != '%':
                comment = False

        n_rows, n_cols, n_cells = line.strip().split()
        n_rows = int(n_rows)
        n_cols = int(n_cols)
        self._n_cells = int(n_cells)

        expected_rows = len(self.rows)
        expected_cols = len(self.columns)
        if n_rows != expected_rows:
            raise Exception(f"Row counts mismatch ({expected_rows}/{n_rows})")
        if n_cols != expected_cols:
            raise Exception(f"Column counts mismatch ({expected_cols}/{n_cols})")

    def next(self):
        """Gets the next line from the table.
        
        This method returns a tuple containing:
        - the row index,
        - the column index,
        - the value.
        """

        if self._handle is None:
            self.open()

        if self._progress_cb:
            line_number = 0
            progress = 0
            one_percent = int(self._n_cells / 100)

        for line in self._handle:
            x, y, v = line.rstrip().split()
            x = self._rows[int(x) - 1]
            y = self._columns[int(y) - 1]
            v = float(v)

            if self._progress_cb:
                line_number += 1
                if line_number % one_percent == 0:
                    progress += 1
                    self._progress_cb(progress)

            yield (x, y, v)

    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        return self.next()
