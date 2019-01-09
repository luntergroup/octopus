#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from math import log, log10, ceil

sns.set(style="darkgrid")

base_palette = {'A': sns.color_palette()[2], 'C': sns.color_palette()[0], 'G': sns.color_palette()[5], 'T': sns.color_palette()[3],}

def estimate_heterozygosity(polymorphisms, opportunities):
    p = polymorphisms / max(2 * opportunities, 1)
    return 2 * p * (1 - p)

def aggregate_polymorphisms(indel_profile_df, aggregators=['platform', 'period', 'periods']):
    if len(aggregators) == 0:
        return indel_profile_df
    polymorphism_attributes = ['platform', 'period', 'periods', 'tract_length', 'motif', 'reference_count', 'reference_span', 'indel_length', 'polymorphisms']
    polymorphism_df = indel_profile_df[polymorphism_attributes]
    aggregation_functions = {'platform': 'first', 'period': 'first', 'periods': 'first', 'reference_count': 'sum', 'reference_span': 'sum', 'polymorphisms': 'sum', 'tract_length': 'first'}
    if 'indel_length' not in aggregators:
        aggregation_functions['reference_count'] = 'first'
        aggregation_functions['reference_span'] = 'first'
        aggregation_functions['motif'] = 'first'
        if 'motif' in aggregators:
            polymorphism_df = polymorphism_df.groupby(aggregators).aggregate(aggregation_functions, axis='columns')
        else:
            polymorphism_df = polymorphism_df.groupby(aggregators + ['motif']).aggregate(aggregation_functions, axis='columns')
        aggregation_functions['reference_count'] = 'sum'
        aggregation_functions['reference_span'] = 'sum'
    else:
        aggregation_functions['indel_length'] = 'first'
    if 'motif' in aggregators:
        aggregation_functions['motif'] = 'first'
    polymorphism_df = polymorphism_df.groupby(aggregators).aggregate(aggregation_functions, axis='columns')
    polymorphism_df['polymorphism_rate'] = polymorphism_df.apply(lambda row: row.polymorphisms / max(row.reference_count, 1), axis=1)
    polymorphism_attributes.append('polymorphism_rate')
    polymorphism_df['heterozygosity'] = polymorphism_df.apply(lambda row: estimate_heterozygosity(row.polymorphisms, row.reference_count), axis=1)
    polymorphism_attributes.append('heterozygosity')
    if 'motif' not in aggregators:
        del(polymorphism_attributes[polymorphism_attributes.index('motif')])
    if 'indel_length' not in aggregators:
        del(polymorphism_attributes[polymorphism_attributes.index('indel_length')])
    return polymorphism_df[polymorphism_attributes]

def aggregate_errors(df, aggregators=['platform', 'period', 'periods']):
    if len(aggregators) == 0:
        return df
    errors_attributes = ['platform', 'period', 'periods', 'tract_length', 'motif', 'indel_length', 'errors', 'reads']
    errors_df = df[errors_attributes]
    aggregation_functions = {'platform': 'first', 'period': 'first', 'periods': 'first', 'errors': 'sum', 'reads': 'sum', 'tract_length': 'first'}
    if 'indel_length' not in aggregators:
        aggregation_functions['motif'] = 'first'
        aggregation_functions['reads'] = 'first'
        if 'motif' in aggregators:
            errors_df = errors_df.groupby(aggregators).aggregate(aggregation_functions, axis='columns')
        else:
            errors_df = errors_df.groupby(aggregators + ['motif']).aggregate(aggregation_functions, axis='columns')
        aggregation_functions['reads'] = 'sum'
    else:
        aggregation_functions['indel_length'] = 'first'
    if 'motif' in aggregators:
        aggregation_functions['motif'] = 'first'
    errors_df = errors_df.groupby(aggregators).aggregate(aggregation_functions, axis='columns')
    errors_df['error_rate'] = errors_df.apply(lambda row: row.errors / max(row.reads, 1), axis=1)
    errors_attributes.append('error_rate')
    if 'motif' not in aggregators:
        del(errors_attributes[errors_attributes.index('motif')])
    if 'indel_length' not in aggregators:
        del(errors_attributes[errors_attributes.index('indel_length')])
    return errors_df[errors_attributes]

def read_indel_profile(csv, platform, model=None):
    result = pd.read_csv(csv)
    result['platform'] = platform
    result['polymorphism_rate'] = result.apply(lambda row: row.polymorphisms / max(row.reference_count, 1), axis=1)
    result['heterozygosity'] = result.apply(lambda row: estimate_heterozygosity(row.polymorphisms, row.reference_count), axis=1)
    result['error_rate'] = result.apply(lambda row: row.errors / max(row.reads, 1), axis=1)
    result['tract_length'] = result.apply(lambda row: row.period * row.periods, axis=1)
    if model is not None:
        result['model'] = model
    return result

def read_indel_profiles(csvs, model=None):
    dfs = []
    for platform, csvs in csvs.items():
        for csv in csvs:
            dfs.append(read_indel_profile(csv, platform, model=model))
    return pd.concat(dfs)

def make_polymorphism_summaries(indel_profile_df):
    result = {}
    result['complex'] = indel_profile_df.query('motif == "N" and polymorphisms > 0')
    result['homopolymer'] = indel_profile_df.query('period == 1 and periods > 1 and motif!="N"')
    result['dinucleotide'] = indel_profile_df.query('period == 2 and periods > 1 and motif not in ["AA", "CC", "GG", "TT"]')
    result['trinucleotide'] = indel_profile_df.query('period == 3 and periods > 0 and motif not in ["AAA", "CCC", "GGG", "TTT"]')
    result['period'] = aggregate_polymorphisms(indel_profile_df)
    result['length'] = aggregate_polymorphisms(indel_profile_df, ['platform', 'periods', 'period', 'indel_length'])
    motif_aggregators = ['platform', 'periods', 'period', 'motif']
    result['motif'] = aggregate_polymorphisms(indel_profile_df.query('motif != "N"'), motif_aggregators)
    result['homopolymer-motif'] = aggregate_polymorphisms(result['homopolymer'], motif_aggregators)
    result['dinucleotide-motif'] = aggregate_polymorphisms(result['dinucleotide'], motif_aggregators)
    result['trinucleotide-motif'] = aggregate_polymorphisms(result['trinucleotide'], motif_aggregators)
    return result

def make_error_summaries(indel_profile_df):
    result = {}
    result['complex'] = indel_profile_df.query('motif == "N" and errors > 0')
    result['homopolymer'] = indel_profile_df.query('period == 1 and periods > 0 and motif != "N"')
    result['dinucleotide'] = indel_profile_df.query('period == 2 and periods > 0 and motif not in ["AA", "CC", "GG", "TT"]')
    result['trinucleotide'] = indel_profile_df.query('period == 3 and periods > 0 and motif not in ["AAA", "CCC", "GGG", "TTT"]')
    result['period'] = aggregate_errors(indel_profile_df)
    result['length'] = aggregate_errors(indel_profile_df, ['platform', 'periods', 'period', 'indel_length'])
    motif_aggregators = ['platform', 'periods', 'period', 'motif']
    result['homopolymer-motif'] = aggregate_errors(result['homopolymer'], motif_aggregators)
    result['dinucleotide-motif'] = aggregate_errors(result['dinucleotide'], motif_aggregators)
    result['trinucleotide-motif'] = aggregate_errors(result['trinucleotide'], motif_aggregators)
    return result

def read_indel_profile_summaries(indel_profile_df):
    result = {}
    result['polymorphism'] = make_polymorphism_summaries(indel_profile_df)
    result['error'] = make_error_summaries(indel_profile_df)
    return result

def read_indel_profile_and_summary(csv, platform, model=None):
    profile_df = read_indel_profile(csv, platform, model=model)
    result = read_indel_profile_summaries(profile_df)
    result['full'] = profile_df
    return result

def read_indel_profile_and_summaries(csvs, model=None):
    profile_df = read_indel_profiles(csvs, model=model)
    result = read_indel_profile_summaries(profile_df)
    result['full'] = profile_df
    return result

def make_octopus_indel_error_model(profile_df):
    # Rules: open[i] <= open[i + 1] + extend[i + 1]
    df = profiles_df['error']['period'].query('platform == "HiSeq-GIAB-HG005" and period == 3 and periods <= 50')
    df['phred_error'] = df.apply(lambda row: -10 * log10(row['error_rate'] + 1e-100), axis=1)
    penalties = [53, 53]
    for periods in range(2, 45):
        penalties.append(float(list(df.query('periods == ' + str(periods))['phred_error'])[0]))
    print(','.join([str(p) for p in [max(int(ceil(penalty)) + 2, 3) for penalty in penalties]]))

def main(options):
    profiles_df = read_indel_profile_and_summaries()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-P','--profile',
                        type=str,
                        required=True,
                        help='Octopus data profile')
    parser.add_argument('-O','--output',
                        type=str,
                        required=True,
                        help='Output prefix')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
