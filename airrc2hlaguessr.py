#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
airrc2hlaguessr.py

Quick script to convert TCR rearrangement files in AIRR-C format into HLAGuessr-compatible files.
"""


import argparse
import sys
import collections as coll
import pandas as pd

__version__ = '0.1.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(description="airrc2hlaguessr v" + str(__version__) + '\n' +
                                                 ": Convert AIRR-C TCR repertoire files to the HLAGuessr input format.")

    parser.add_argument('-ia', '--ignore_ambiguous', action='store_true', required=False, default=False,
                        help="Optional flag to ignore TCRs with ambiguous V calls.\n"
                             "(Default behaviour is to include their CDR3s with each gene family redundantly.)")

    parser.add_argument('-z', '--compress', action='store_true', required=False, default=False,
                        help="Optional flag to gzip compress the output.")

    parser.add_argument('-kd', '--keep_duplicates', action='store_true', required=False, default=False,
                        help="Optional flag to keep duplicate TCRs. Useful when processing very large datasets. ")

    parser.add_argument('-ts', '--truncate_str', type=str, required=False,
                        help='Optional string to truncate filenames for in output.'
                             '\nUseful for ensuring same name used for both alpha and beta files.')

    parser.add_argument('-in', '--in_files', type=str, required=False, default='',
                        help="Optional field to input paths to files manually. Comma delimited.")

    parser.add_argument('-c', '--chain_filter', type=str, required=False, default='',
                        help="Optional field to filter TCR chains. Only 'A' or 'B' allowed. ")

    parser.add_argument('-n', '--name_override', type=str, required=False, default='',
                        help="Optional field to override the name used in the ID field,"
                             " instead of inferring from filename (overrides -ts).")

    return parser.parse_args()


def list_to_df(input_list, headers, rename):
    """
    Convert a list to a (long) dataframe. Note that first entry becomes the index if chosen
    :param input_list: List of list entries (with each position in each list corresponding to a column)
    :param headers: List of column headers. First column should be unique, becoming the rownames, if rename = True
    :param rename: Option to rename row IDs by first colum
    :return: sorted pandas dataframe
    """
    df = pd.DataFrame(input_list)
    df = df.rename(index=str, columns=dict(list(zip(list(range(len(headers))), headers))))
    if rename == True:
        df = df.set_index(headers[0], drop=True)
    return df


if __name__ == '__main__':

    input_args = vars(args())

    # If file paths provided explicitly, use them; otherwise take piped input
    if input_args['in_files']:
        in_file_list = input_args['in_files'].split(',')
    else:
        in_file_list = sys.stdin.read().rsplit()

    if not in_file_list:
        raise IOError("No input files detected. Please pipe a valid file list, or use the -in field to specify. ")

    if input_args['chain_filter']:
        if input_args['chain_filter'][-1].upper() == 'A':
            keep_chains = ['A', 'D']
        elif input_args['chain_filter'][-1].upper() == 'B':
            keep_chains = ['B']
        else:
            raise IOError("Inappropriate chain filter selected: only 'TRA'/'TRB', or 'A'/'B' are valid. ")
    else:
        keep_chains = ''

    for pf in in_file_list:
        print("Processing " + pf + "...")
        try:
            df = pd.read_csv(pf, compression='infer', sep='\t')
            nam = pf[:pf.index('.tsv')].split('/')[-1]
            line_nam = nam
            if input_args['truncate_str']:
                truncate_index = nam.find(input_args['truncate_str'])
                if truncate_index > 0:
                    line_nam = line_nam[:truncate_index]
                else:
                    print("Truncate string provided (" + input_args['truncate_str'] +
                          ") which doesn't appear in filename - ignoring...")
            out_nam = nam + '_hlaguessr.tsv'

            if input_args['name_override']:
                line_nam = input_args['name_override']

            if input_args['compress']:
                out_nam += '.gz'

            out_dat = []
            tcrs = coll.Counter()
            for tcr in df.index:
                row_dat = df.loc[tcr]

                # Skip non-productives
                if row_dat['productive'] == 'F' or not isinstance(row_dat['junction_aa'], str):
                    continue
                elif len(row_dat['junction_aa']) == 0:
                    continue

                v_fam = list(set([x.split('*')[0].split('-')[0] for x in row_dat['v_call'].split(',')]))

                # By default, only output rearrangements with unambiguous V calls
                if len(v_fam) == 1:
                    out_dat.append([row_dat['junction_aa'], v_fam[0], line_nam])

                # If the appropriate input flag is set, write out a line each for ambiguous V calls
                elif len(v_fam) > 1 and not input_args['ignore_ambiguous']:
                    for specific_v in v_fam:
                        out_dat.append([row_dat['junction_aa'], specific_v, line_nam])

                else:
                    continue

            out_dat = list_to_df(out_dat, ['cdr3aa', 'v_family', 'Patient'], False)
            if not input_args['keep_duplicates']:
                out_dat = out_dat.drop_duplicates().sort_values(by=['cdr3aa', 'v_family'], axis=0)
            if keep_chains:
                out_dat = out_dat.loc[out_dat['v_family'].str[2].isin(keep_chains)]
            out_dat.to_csv(out_nam, sep='\t', index=False, compression='infer')

        except Exception:
            print("Failed to read in file " + pf + " - skipping...")
