# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

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
