# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import os.path
import subprocess
from shutil import copyfile

import click
from click_shell.core import make_click_shell
from pandas import concat, read_csv, DataFrame


grep_svm_script = """
SVM_DIR=$1
REGEX=$2

find $SVM_DIR -type f | while read FILE ; do
    BASENAME=${FILE##*/}
    case "$BASENAME" in
    FBrf*)
        echo -e "$BASENAME\\t$FILE" >> $$.fulltext-by-fbrf
        ;;
    *)
        echo -e "$BASENAME\\t$FILE" >> $$.fulltext-by-pmid
    esac
done

for REF in $REFS ; do
    FBRF=${REF%:*}
    PMID=${REF#*:}
    FULLTEXT_PATH=$(cat $$.fulltext-by-fbrf | grep ^$FBRF | cut -f2)
    if [ -z "$FULLTEXT_PATH" ]; then
        FULLTEXT_PATH=$(cat $$.fulltext-by-pmid | grep ^$PMID | cut -f2)
    fi

    if [ -z "$FULLTEXT_PATH" ]; then
        echo -e "$FBRF\\tFull-text not available"
    elif [ ! -f "$FULLTEXT_PATH" ]; then
        echo -e "$FBRF\\tFull-text not available"
    else
        NMATCH=$(grep -i -c -E "$REGEX" $FULLTEXT_PATH)
        echo -e "$FBRF\\t$NMATCH"
    fi
done

rm $$.fulltext-by-{fbrf,pmid}
"""


class DiscoverContext(object):
    """A helper object for all discovery operations."""

    def __init__(self, config, database):
        self._config = config
        self._database = database

    @property
    def flybase_table_file(self):
        return self._config.get(
            'discovery', 'flybase_papers_list', fallback='flybase-papers.tsv'
        )

    @property
    def scea_dataset_file(self):
        return self._config.get('discovery', 'scea_dataset_list')

    @property
    def scea_staging_dataset_file(self):
        return self._config.get('discovery', 'scea_staging_list')

    @property
    def cursor(self):
        return self._database.cursor

    def filter_out_known_datasets(self, pmids):
        table = read_csv(self.scea_dataset_file, dtype=str)
        results = [p for p in pmids if p not in table['PubMed ID'].values]

        table = read_csv(self.scea_staging_dataset_file, dtype=str)
        return [p for p in results if p not in table['PMID'].values]

    def get_dataset_fbrfs(self):
        """Get all references that have been flagged with 'dataset'."""

        query = f'''SELECT pub.uniquename
                    FROM
                             pub
                        JOIN pubprop USING (pub_id)
                        JOIN cvterm  ON pubprop.type_id = cvterm.cvterm_id
                        JOIN cv      USING (cv_id)
                    WHERE
                             cv.name         LIKE 'pubprop type'
                         AND cvterm.name     LIKE 'harv_flag'
                         AND pubprop.value   LIKE 'dataset'
                         AND pub.is_obsolete = 'f';'''
        self.cursor.execute(query)
        return self.cursor.fetchall()

    def get_internalnotes_for_fbrf(self, fbrf):
        """Get the internal notes for a given reference."""

        query = f'''SELECT pubprop.value
                    FROM
                             pub
                        JOIN pubprop USING (pub_id)
                        JOIN cvterm  ON pubprop.type_id = cvterm.cvterm_id
                        JOIN cv      USING (cv_id)
                    WHERE
                             cv.name        LIKE 'pubprop type'
                         AND cvterm.name    LIKE 'internalnotes'
                         AND pub.uniquename LIKE '{fbrf}';'''
        self.cursor.execute(query)
        return self.cursor.fetchall()

    def get_dataset_accessions_for_fbrf(self, fbrf):
        """Get dataset accession numbers associated with a reference."""

        notes = self.get_internalnotes_for_fbrf(fbrf)
        accessions = []
        for note in notes:
            for line in note[0].split('\n'):
                if not line.startswith('Dataset:'):
                    continue
                line = line[8:].strip().split('.')[0]
                for item in line.split(','):
                    item = item.strip()
                    if item.startswith('and '):
                        item = item[3:].strip()
                    if item == 'pheno':
                        continue
                    accessions.append(item)
        return ','.join(accessions)

    def get_pmid_for_fbrf(self, fbrf):
        """Get the PubMed ID for a given reference."""

        query = f'''SELECT dbxref.accession
                    FROM
                             pub
                        JOIN pub_dbxref USING (pub_id)
                        JOIN dbxref     USING (dbxref_id)
                        JOIN db         USING (db_id)
                    WHERE
                             db.name        LIKE 'pubmed'
                         AND pub.uniquename LIKE '{fbrf}';'''
        self.cursor.execute(query)
        return ','.join([a for a, in self.cursor.fetchall()])

    def get_citation_for_fbrf(self, fbrf):
        """Get a text citation for a given reference."""

        query = f'''SELECT pubauthor.surname, pubauthor.rank, pub.pyear
                    FROM
                             pubauthor
                        JOIN pub USING (pub_id)
                    WHERE
                             pub.uniquename LIKE '{fbrf}'
                    ORDER BY pubauthor.rank;'''
        self.cursor.execute(query)
        res = self.cursor.fetchall()
        year = res[0][2]
        n = len(res)
        if n == 1:
            return f'{res[0][0]}, {year}'
        elif n == 2:
            return f'{res[0][0]} and {res[1][0]}, {year}'
        else:
            return f'{res[0][0]} et al., {year}'

    def grep_fulltexts(self, references):
        """Search the full-texts archive for references to scRNAseq.

        :param references: a list of tuple (FBrf, PMID)
        :return: a dictionary whose keys are FBrfs and values are numbers
                 of references to scRNAseq
        """

        hostname = self._config.get('textmining', 'host')
        directory = self._config.get('textmining', 'directory')
        pattern = self._config.get('textmining', 'pattern')

        script = (
            'REFS="'
            + ' '.join([f'{fbrf}:{pmid}' for fbrf, pmid in references])
            + '"\n'
            + grep_svm_script
        )
        command = f"bash -s {directory} '{pattern}'"
        ssh = subprocess.Popen(
            ['ssh', hostname, command],
            shell=False,
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        output, _ = ssh.communicate(script)

        res = {}
        for line in output.split('\n'):
            line = line.strip()
            if '\t' not in line:
                continue
            fbrf, nmatches = line.split('\t')
            if nmatches == 'Full-text not available':
                res[fbrf] = -1
            else:
                res[fbrf] = int(nmatches)
        return res


@click.group(invoke_without_command=True)
@click.pass_context
def discover(ctx):
    """Access the dataset discovery commands."""

    ctx.obj = DiscoverContext(ctx.obj.config, ctx.obj.database)

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-discover> ")
        shell.cmdloop()


def findnew_impl(obj, table, fbrfs):
    """Update the given table with data for the specified references.

    This is the core of the findnew/checknew commands.
    """

    if table is None:
        table = DataFrame(
            columns=[
                'FBrf',
                'PMID',
                'Citation',
                'Accessions',
                'Mentions',
                'Confirmed',
                'Organ/tissue',
                'Comments',
            ],
            dtype=str,
        )

    click.echo("Fetching additional data...")
    newrows = []
    with click.progressbar(fbrfs) as bar:
        for fbrf in bar:
            if fbrf in table['FBrf'].values:
                continue

            pmid = obj.get_pmid_for_fbrf(fbrf)
            accessions = obj.get_dataset_accessions_for_fbrf(fbrf)
            citation = obj.get_citation_for_fbrf(fbrf)
            newrows.append(
                {
                    'FBrf': fbrf,
                    'PMID': pmid,
                    'Citation': citation,
                    'Accessions': accessions,
                }
            )

    click.echo(f"New dataset-flagged references: {len(newrows)}")
    if len(newrows) > 0:
        table = concat([table, DataFrame(data=newrows)])

        click.echo("Querying the fulltext archive...")
        queries = table.loc[table['Mentions'].isna(), ['FBrf', 'PMID']].values
        mentions = obj.grep_fulltexts(queries)
        pos = 0
        nas = 0
        for fbrf, nmatch in mentions.items():
            if nmatch > 0:
                pos += 1
            elif nmatch == -1:
                nas += 1
                nmatch = 'Full-text not available'
            table.loc[table['FBrf'] == fbrf, 'Mentions'] = nmatch
        click.echo(f"References matching the scRNAseq pattern: {pos}")
        click.echo(f"References without full text: {nas}")

    return len(newrows), table


@discover.command()
@click.option(
    '--filename',
    '-f',
    default=None,
    metavar='FILENAME',
    help="""Use the given table file instead of the one
                      specified in the configuration.""",
)
@click.option(
    '--output',
    '-o',
    default=None,
    metavar='FILENAME',
    help="""Write the updated table to the specified file.
            The default is to write to the original file.""",
)
@click.pass_obj
def findnew(obj, filename, output):
    """Find new FlyBase references that may be about scRNAseq."""

    if filename is None:
        filename = obj.flybase_table_file

    table = None
    if os.path.exists(filename):
        table = read_csv(filename, sep='\t', dtype=str)

    click.echo("Fetching dataset-flagged references...")
    fbrfs = [f[0] for f in obj.get_dataset_fbrfs()]

    n, table = findnew_impl(obj, table, fbrfs)
    if n > 0:
        if output is None:
            output = filename
            copyfile(output, f'{output}.bak')
        table.to_csv(output, index=False, sep='\t')


@discover.command()
@click.option(
    '--output',
    '-o',
    type=click.File('w'),
    default='-',
    metavar='FILENAME',
    help="""Write the result to the specified file.
            The default is to write to standard output."""
)
@click.argument('fbrfs', nargs=-1)
@click.pass_obj
def checknew(obj, output, fbrfs):
    """Check whether the provided references may be about scRNAseq."""

    n, table = findnew_impl(obj, None, fbrfs)
    if n > 0:
        table.to_csv(output, index=False, sep='\t')


@discover.command()
@click.option(
    '--output',
    '-o',
    type=click.File('w'),
    default='-',
    help="Send output to the specified file.",
)
@click.pass_obj
def toscea(obj, output):
    """Make a list of datasets that are new to the SCEA."""

    table = read_csv(obj.flybase_table_file, sep='\t', dtype=str)
    subset = table[table['Confirmed'] == 'yes']
    unknown_pmids = obj.filter_out_known_datasets(subset['PMID'].values)

    subset = table[table['PMID'].isin(unknown_pmids)]
    subset = table.loc[
        table['PMID'].isin(unknown_pmids),
        ['FBrf', 'PMID', 'Accessions', 'Citation', 'Organ/tissue', 'Comments'],
    ]
    subset.to_csv(output, sep='\t', index=False)
