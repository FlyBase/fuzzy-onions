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

from psycopg2 import connect


class DatabaseHelper(object):
    """A helper object to connect to and use the FlyBase database."""

    def __init__(self, config):
        self._config = config
        self._conn = None
        self._cursor = None

    @property
    def connection(self):
        if self._conn is None:
            host = self._config.get('chado', 'host', fallback='chado.flybase.org')
            name = self._config.get('chado', 'database', fallback='flybase')
            user = self._config.get('chado', 'user', fallback='flybase')
            pswd = self._config.get('chado', 'password', fallback=None)
            name = self._config.get('chado', 'database', fallback='latest')
            if name == 'latest':
                name = self._get_latest_database(host, user, pswd)
            self._conn = connect(host=host, database=name, user=user, password=pswd)
        return self._conn

    @property
    def cursor(self):
        if self._cursor is None:
            self._cursor = self.connection.cursor()
        return self._cursor

    def close(self):
        if self._cursor is not None:
            self._cursor.close()
        if self._conn is not None:
            self._conn.close()

    def _get_latest_database(self, host, user, password):
        query = f'''SELECT datname
                    FROM
                             pg_database
                    WHERE
                             datistemplate = false
                        AND  datname LIKE 'fb_20%'
                    ORDER BY datname DESC
                    LIMIT 1;'''
        with connect(
            host=host, database='postgres', user=user, password=password
        ) as tmp:
            with tmp.cursor() as cursor:
                cursor.execute(query)
                return cursor.fetchone()[0]
