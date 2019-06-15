#!/usr/bin/env python3

import argparse
import pandas as pd
from math import log10, ceil

def aggregate_errors(df, aggregators=['period', 'periods']):
    if len(aggregators) == 0:
        return df
    if len(df) == 0:
        return df
    errors_attributes = ['period', 'periods', 'tract_length', 'motif', 'indel_length', 'errors', 'reads']
    errors_df = df[errors_attributes]
    aggregation_functions = {'period': 'first', 'periods': 'first', 'errors': 'sum', 'reads': 'sum', 'tract_length': 'first'}
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
    errors_df.reset_index(drop=True, inplace=True)
    errors_df = errors_df.groupby(aggregators).aggregate(aggregation_functions, axis='columns')
    errors_df['error_rate'] = errors_df.apply(lambda row: row.errors / max(row.reads, 1), axis=1)
    errors_attributes.append('error_rate')
    if 'motif' not in aggregators:
        del(errors_attributes[errors_attributes.index('motif')])
    if 'indel_length' not in aggregators:
        del(errors_attributes[errors_attributes.index('indel_length')])
    return errors_df[errors_attributes]

def read_indel_profile(csv):
    result = pd.read_csv(csv)
    result['error_rate'] = result.apply(lambda row: row.errors / max(row.reads, 1), axis=1)
    result['tract_length'] = result.apply(lambda row: row.period * row.periods, axis=1)
    return result

def make_error_summaries(indel_profile_df):
    result = {}
    result['complex'] = indel_profile_df.query('motif == "N" and errors > 0')
    result['homopolymer'] = indel_profile_df.query('period == 1 and periods > 0 and motif != "N"')
    result['dinucleotide'] = indel_profile_df.query(
        'period == 2 and periods > 0 and motif not in ["AA", "CC", "GG", "TT"]')
    result['trinucleotide'] = indel_profile_df.query(
        'period == 3 and periods > 0 and motif not in ["AAA", "CCC", "GGG", "TTT"]')
    result['tetranucleotide'] = indel_profile_df.query(
        'period == 4 and periods > 0 and motif not in ["AAAA", "ACAC", "AGAG", "ATAT", "CACA", "CCCC", "CGCG", "CTCT", "GAGA", "GCGC", "GGGG", "GTGT", "TATA", "TCTC", "TGTG", "TTTT"]')
    result['pentanucleotide'] = indel_profile_df.query(
        'period == 5 and periods > 0 and motif not in ["AAAAA", "CCCCC", "GGGGG", "TTTTT"]')
    result['period'] = aggregate_errors(indel_profile_df)
    result['length'] = aggregate_errors(indel_profile_df, ['periods', 'period', 'indel_length'])
    motif_aggregators = ['periods', 'period', 'motif']
    result['homopolymer-motif'] = aggregate_errors(result['homopolymer'], motif_aggregators)
    result['dinucleotide-motif'] = aggregate_errors(result['dinucleotide'], motif_aggregators)
    result['trinucleotide-motif'] = aggregate_errors(result['trinucleotide'], motif_aggregators)
    result['tetranucleotide-motif'] = aggregate_errors(result['tetranucleotide'], motif_aggregators)
    result['pentanucleotide-motif'] = aggregate_errors(result['pentanucleotide'], motif_aggregators)
    return result

def read_indel_error_profile_summaries(indel_profile_df):
    return make_error_summaries(indel_profile_df)

def read_indel_error_profile_and_summary(csv, model=None):
    return read_indel_error_profile_summaries(read_indel_profile(csv))

def rate_to_phred(rate):
    return -10 * log10(rate + 1e-100)

def complex_indel_penalty(profile_df):
    return rate_to_phred(float(list(profile_df['period'].query('period == 0')['error_rate'])[0]))

def get_repeat_error_df(profile_df, pattern, max_periods):
    profile_index = 'period'
    query_condition = "periods <= " + str(max_periods)
    if type(pattern) == int:
        query_condition += " and period == " + str(pattern)
    else:
        if len(pattern) == 1:
            motif = pattern[0]
            motif_len = len(motif)
            query_condition += " and motif == '" + motif + "'"
        else:
            motif_len = len(pattern[0])
            query_condition += " and (motif == '" + pattern[0] + "'"
            for motif in pattern[1:]:
                " or motif == '" + motif + "'"
            query_condition += ")"
        if motif_len == 1:
            profile_index = 'homopolymer-motif'
        elif motif_len == 2:
            profile_index = 'dinucleotide-motif'
        elif motif_len == 3:
            profile_index = 'trinucleotide-motif'
    result = profile_df[profile_index].query(query_condition).copy()
    result['phred_error'] = result['error_rate'].apply(rate_to_phred)
    return result

def make_empirical_indel_error_model_helper(profile_df, pattern, max_periods=50):
    repeat_error_df = get_repeat_error_df(profile_df, pattern, max_periods)
    complex_penalty = complex_indel_penalty(profile_df)
    penalties = [complex_penalty, complex_penalty]
    for periods in range(2, max_periods):
        try:
            penalties.append(float(list(repeat_error_df.query('periods == ' + str(periods))['phred_error'])[0]))
        except IndexError:
            break
    return [max(int(ceil(penalty)) + 2, 3) for penalty in penalties]

def max_lt(seq, val):
    return max(v for v in seq if v < val)

def smooth_empirical_model(open_model, extend_model=10):
    if len(open_model) < 3:
        return open_model
    # Rules: open[i] <= open[i + 1] + extend[i + 1]
    head_penalties = open_model[:2]
    max_penalty = max_lt(open_model[2:], 50)
    result = [max_penalty, max_penalty]
    for i, penalty in enumerate(open_model[2:], 2):
        if i < len(open_model) - 1 and penalty < open_model[i + 1]:
            penalty = int((result[i - 1] + min(result[i - 1], open_model[i + 1])) / 2)
        result.append(max(min(penalty, result[i - 1]), result[i - 1] - extend_model))
    result[:2] = head_penalties
    return result

def make_octopus_indel_error_model(profile_df):
    result = {}
    for pattern in [('A', 'T'), ('C', 'G'), 2, 3, 4, 5]:
        model = smooth_empirical_model(make_empirical_indel_error_model_helper(profile_df, pattern))
        if type(pattern) == int:
            result[pattern * 'N'] = model
        else:
            for motif in pattern:
                result[motif] = model
    if 'N' not in result:
        result['N'] = result['A']
    return result

def write_error_model(error_model, out_path):
    with open(out_path, 'w') as file:
        for motif, model in error_model.items():
            file.write(motif + ":" + ','.join([str(p) for p in model]) + '\n')

def main(options):
    profiles_df = read_indel_error_profile_and_summary(options.profile)
    error_model = make_octopus_indel_error_model(profiles_df)
    write_error_model(error_model, options.output)

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
