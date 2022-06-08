from argparse import ArgumentParser
import os
import json
import pandas as pd
import pprint
from functools import reduce
import numpy as np

def createArgParser():
    # Create argument parser
    description = "Parse HLA-HD final results files and determine consensus HLA type for each HLA gene."
    argparser = ArgumentParser(description = description)
    argparser.add_argument('-s1', '--sample_1_results', required=True, dest='sample_1_results', type=str, 
                           help='Path to sample 1 HLA-HD final results file')
    argparser.add_argument('-s2', '--sample_2_results', required=True, dest='sample_2_results', type=str, 
                           help='Path to sample 2 HLA-HD final results file')
    argparser.add_argument('-s3', '--sample_3_results', required=False, dest='sample_3_results', type=str, 
                           help='Path to (optional) sample 3 HLA-HD final results file')
    argparser.add_argument('patient_id', type=str, 
                           help='Patient ID')
    
    return argparser

def parseHLAHD(results_file):
    # Check that the file path is valid
    if not os.path.exists(results_file):
        raise Exception("HLA-HD results file does not exist")

    # Open the results file
    with open(results_file, "r") as f:
        allele_dict_short = {}
        allele_dict_full = {}

        for line in f:
            line = line.strip().split('\t')

            # Gene name
            gene = line[0]
            # All listed alleles (can be more than 2)
            alleles = line[1:]
            alleles_edited = []

            for i, a in enumerate(alleles):
                if a not in ['Not typed', '-']:
                    # Strip the gene name from the front
                    num = a.replace('HLA-' + gene + '*', '')
                    # Separate digit accuracies into list for each allele
                    alleles_edited.append(num)
                elif a == '-':
                    # Homozygous for this allele, so copy the previous allele
                    alleles_edited.append(alleles_edited[i-1])
                else:
                    alleles_edited.append(a)

            allele_dict_full[gene] = alleles_edited

    return allele_dict_full

def df_to_dict(consensus_df):
    consensus_dict = {}

    for i, row in consensus_df.iterrows():
        consensus_dict['HLA-' + row['hla_gene']] = row[1]

    return consensus_dict

def clinical_hla_json(hla_dict, filename):
    clinical_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1']
    clinical = {k: hla_dict[k] for k in clinical_genes}

    # for k in clinical.keys():
    #     truncated = []
    #     for a in clinical[k]:
    #         a_trunc = ':'.join(a.split(':')[:2])
    #         truncated.append(a_trunc)
    #     clinical[k] = truncated




    write_to_json(clinical, filename)

    return clinical

def write_to_json(allele_dict, filename):
    with open(filename, "w") as out: 
        json.dump(allele_dict, out, indent=4)
        out.close()
    return

def format_df(consensus_df):
    consensus_df['hla_gene'] = 'HLA-' + consensus_df['hla_gene']

    formatted_df = pd.DataFrame(columns=['Allele 1', 'Allele 2'])
    formatted_df.index.name = 'HLA Gene'

    for i, row in consensus_df.iterrows():
        for allele in row[1]:
            formatted_df.loc[row['hla_gene'], 'Allele 1'] = row[1][0]
            formatted_df.loc[row['hla_gene'], 'Allele 2'] = row[1][1]

    return formatted_df

def trunc_hla_txt(hla_dict, filename):
    allele_list = []

    for gene in hla_dict.keys():
        if gene not in ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G']:
            pvac_gene = gene.replace('HLA-', '')
        else:
            pvac_gene = gene

        for allele in hla_dict[gene]:
            if allele in ['Not typed', 'No consensus']:
                continue
            a_trunc = ':'.join(allele.split(':')[:2])
            full_allele = pvac_gene + '*' + a_trunc

            if not full_allele in allele_list:
                allele_list.append(full_allele)

    allele_string = '\'' + '\',\''.join(allele_list) + '\''
    print(allele_string)

    write_to_string(allele_string, filename)

    return allele_string

def write_to_string(allele_string, filename):
    with open(filename, "w") as out:
        out.write(allele_string)
        out.close
    return

def dict_to_df(sample_dict):
    df = pd.DataFrame([sample_dict]).transpose()
    df.index.name = "hla_gene"
    df.columns = ['alleles']
    return df

def find_mode(df, col):
    counts = df[col].value_counts()
    mode_values = counts[counts.eq(counts.max())]

    if len(mode_values) == 1:
        # There was a clear winner (no tie)
        return mode_values.index[0]
    elif len(mode_values) > 1:
        for i, m in enumerate(mode_values):
            if m >= 2:
                # Sufficient evidence since found in at least 2 samples
                return mode_values.index[i]

        # Insufficient evidence since only 1 sample has this value
        return None

def truncate_alleles(alleles):
    trunc_alleles = []

    for a in alleles:
        if a == "Not typed":
            trunc_alleles.append(a)
        else:
            a_split = a.split(":")
            trunc_alleles.append(":".join(a_split[:2]))

    return trunc_alleles

def cross(gene, gene_df, consensus_df):
    # Get number of samples (2 or 3)
    n_samples = len(gene_df.index)

    # Find index of this gene in the consensus df
    i = consensus_df.index[consensus_df['hla_gene'] == gene][0]

    # Create a copy of the gene_df to track which alleles have been used and which remain
    # Note that the df has to be copied this way because of list (mutable object) in pandas df
    remaining_df = pd.DataFrame(columns=['alleles'])
    if n_samples == 3:
        for n in [0, 1, 2]:
            remaining_df = remaining_df.append({'alleles' : list(gene_df['alleles'])[n]}, ignore_index=True)
    elif n_samples == 2:
        for n in [0, 1]:
            remaining_df = remaining_df.append({'alleles' : list(gene_df['alleles'])[n]}, ignore_index=True)

    if n_samples == 3:
        # Cross check each allele called from the first sample
        for j, a in enumerate(list(gene_df['alleles'])[0]):
            if a in list(remaining_df['alleles'])[1] or a in list(remaining_df['alleles'])[2]:
                # Sufficient evidence so add to consensus df
                consensus_df['consensus_alleles'].iloc[i].append(a)
                # Remove the used allele
                m = list(remaining_df['alleles'])[0].index(a)
                remaining_df['alleles'][0] = [al for i, al in enumerate(list(remaining_df['alleles'])[0]) if i != m]

                if a in list(remaining_df['alleles'])[1]:
                    # Get index of first occurence of this allele
                    k = list(remaining_df['alleles'])[1].index(a)
                    # Remove the used allele
                    remaining_df['alleles'][1] = [al for i, al in enumerate(list(remaining_df['alleles'])[1]) if i != k]

                if a in list(remaining_df['alleles'])[2]:
                    # Get index of first occurence of this allele
                    k = list(remaining_df['alleles'])[2].index(a)
                    # Remove the used allele
                    remaining_df['alleles'][2] = [al for i, al in enumerate(list(remaining_df['alleles'])[2]) if i != k]

        for j, a in enumerate(list(remaining_df['alleles'])[1]):
            if a in list(remaining_df['alleles'])[0] or a in list(remaining_df['alleles'])[2]:
                # Sufficient evidence so add to consensus df
                consensus_df['consensus_alleles'].iloc[i].append(a)
                # Remove the used allele
                remaining_df['alleles'][1] = [al for i, al in enumerate(list(remaining_df['alleles'])[1]) if i != j]

                if a in list(remaining_df['alleles'])[0]:
                    # Get index of first occurence of this allele
                    k = list(remaining_df['alleles'])[0].index(a)
                    # Remove the used allele
                    remaining_df['alleles'][0] = [al for i, al in enumerate(list(remaining_df['alleles'])[0]) if i != k]

                if a in list(remaining_df['alleles'])[2]:
                    # Get index of first occurence of this allele
                    k = list(remaining_df['alleles'])[2].index(a)
                    # Remove the used allele
                    remaining_df['alleles'][2] = [al for i, al in enumerate(list(remaining_df['alleles'])[2]) if i != k]

    # Two sample
    elif n_samples == 2:
        for j, a in enumerate(list(gene_df['alleles'])[0]):
            if a in list(remaining_df['alleles'])[1]:
                # Sufficient evidence so add to consensus df
                consensus_df['consensus_alleles'].iloc[i].append(a)
                k = list(remaining_df['alleles'])[1].index(a)
                # Remove used alleles
                remaining_df['alleles'][0] = [al for i, al in enumerate(list(remaining_df['alleles'])[0]) if i != j]
                remaining_df['alleles'][1] = [al for i, al in enumerate(list(remaining_df['alleles'])[1]) if i != k]    


    return consensus_df, remaining_df



def consensus(df_list):
    # Concatenate sample dfs to one long df (multiple rows per gene, 1 for each sample)
    df = pd.concat(df_list)
    # Make 'hla_gene' a columns (rather than index)
    df = df.reset_index()

    # First: look for matching pairs
    df['consensus_alleles'] = df['hla_gene'].map(df.groupby('hla_gene').apply(find_mode, col='alleles'))
    consensus_df = df[['hla_gene', 'consensus_alleles']]
    consensus_df = consensus_df.drop_duplicates(subset=['hla_gene'])

    df = df.sort_values(by=['hla_gene'])

    for i, gene in consensus_df['hla_gene'].iteritems():
        # Find genes for which there was no simple consensus
        if not consensus_df['consensus_alleles'].iloc[i]:
            consensus_df['consensus_alleles'].iloc[i] = []

        gene_df = df[df['hla_gene'] == gene]

        if len(consensus_df['consensus_alleles'].iloc[i]) < 2:
            consensus_df, remaining_df = cross(gene, gene_df, consensus_df)

        # Check 2 field accuracies
        if len(consensus_df['consensus_alleles'].iloc[i]) < 2:  
            # Truncate to 2 field accuracies
            trunc_remaining_df = pd.DataFrame()
            trunc_remaining_df['alleles'] = remaining_df['alleles'].map(truncate_alleles)

            two_field_consensus = pd.DataFrame(columns=['consensus_alleles'])
            # Use the find_mode function to find matching pairs of two-field alleles
            two_field_consensus = find_mode(trunc_remaining_df, 'alleles')

            # Add pair of alleles to consensus df (if found)
            if two_field_consensus:
                consensus_df['consensus_alleles'].iloc[i].extend(two_field_consensus)
            else:
                # Look for single matches
                consensus_df, remaining_df = cross(gene, trunc_remaining_df, consensus_df)

        # Fill in with 'not typed'
        while len(consensus_df['consensus_alleles'].iloc[i]) < 2:
            consensus_df['consensus_alleles'].iloc[i].append('No consensus')

        if len(consensus_df['consensus_alleles'].iloc[i]) > 2:
            # If too many consensus alleles were found (e.g. 3 alleles each of which appeared once in 2 samples) then set as no consensus
            consensus_df['consensus_alleles'].iloc[i] = ['No consensus', 'No consensus']

    return consensus_df

def consensus_matrix(consensus_df, sample_df):
    sample_df = sample_df.reset_index()
    # consensus_df = consensus_df.set_index('hla_gene')
    consensus_df_form = format_df(consensus_df)
    sample_df_form = format_df(sample_df)

    pprint.pprint(sample_df_form)
    pprint.pprint(consensus_df_form)

    print(sample_df_form == consensus_df_form)
    return

# def sample_stats(df):
#     not_typed = 0
#     two_field = 0
#     three_field = 0

#     for i, row in df.iterrows():
#         for a in row['alleles']:
#             if a == 'Not typed':
#                 not_typed = not_typed + 1

#             a_split = a.split(':')
#             if len(a_split) == 2:
#                 two_field = two_field + 1
#             elif len(a_split) == 3:
#                 three_field = three_field + 1

#     return not_typed, two_field, three_field

# def consensus_stats():

# def generate_stats(sample_df_list, consensus_df):
#     for df in sample_df_list:
#         sample_stats(df)
#     return



if __name__ == '__main__':
    argparser = createArgParser()
    args = argparser.parse_args()

    s1_dict = parseHLAHD(args.sample_1_results)
    s1_df = dict_to_df(s1_dict)

    s2_dict = parseHLAHD(args.sample_2_results)
    s2_df = dict_to_df(s2_dict)

    if args.sample_3_results:
        s3_dict = parseHLAHD(args.sample_3_results)
        s3_df = dict_to_df(s3_dict)

        out_suffix = "threeSample"
        df_list = [s1_df, s2_df, s3_df]
    else:
        out_suffix = "twoSample"
        df_list = [s1_df, s2_df]

    # Run consensus algorithm and convert to dictionary
    consensus_df = consensus(df_list)
    
    consensus_dict = df_to_dict(consensus_df)

    clin_filename = args.patient_id + '_' + out_suffix + '_hla.consensus.clinSig.json'
    full_filename = args.patient_id + '_' + out_suffix + '_hla.consensus.json'
    full_txt_filename = args.patient_id + '_' + out_suffix + '_hla.consensus.trunc.txt'
    clin_txt_filename = args.patient_id + '_' + out_suffix + '_hla.consensus.clinSig.trunc.txt'

    # Write json files
    clin_dict = clinical_hla_json(consensus_dict, clin_filename)
    write_to_json(consensus_dict, full_filename)

    # Write text file with typed alleles (for pVACseq input)
    # clinical_hla_txt(clin_dict, txt_filename)
    trunc_hla_txt(consensus_dict, full_txt_filename)
    trunc_hla_txt(clin_dict, clin_txt_filename)

    # cons_matrix = consensus_matrix(consensus_df, s1_df)

    # Write all HLA samples to json files
    s1_filename = args.patient_id + '_' + 'sample1_hla.json'
    s2_filename = args.patient_id + '_' + 'sample2_hla.json'
    write_to_json(df_to_dict(s1_df.reset_index()), s1_filename)
    write_to_json(df_to_dict(s2_df.reset_index()), s2_filename)

    if args.sample_3_results:
        s3_filename = args.patient_id + '_' + 'sample3_hla.json'
        write_to_json(df_to_dict(s3_df.reset_index()), s3_filename)
