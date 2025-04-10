#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# IMPORTS AND FILE LOADING
import pandas as pd
import numpy as np
from tqdm import tqdm
from Bio.Seq import Seq
from Bio import SeqIO 
from Bio.Restriction import *
import re
import argparse
import os
import sys
import importlib.resources as pkg_resources

def parse_and_store_fasta(file_path, v_region_type):
    """
    Parses FASTA file containing ground truth TRAV/TRBV sequences and stores them in a dictionary.
    params
        file_path: Path to the FASTA file.
        v_region_type: Type of V region (TRAV or TRBV).
    return
        Dictionary with sequence names as keys and sequences as values.
    """
    sequence_dict = {}
    with open(file_path, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            name = record.id
            name = name[name.index(v_region_type) : name.index('*')]
            if name not in sequence_dict.keys(): #exclude *02 and higher
                sequence_dict[name] = str(record.seq)
    return sequence_dict

"""
This section loads the reference files needed for the assembly validation.
"""
ref_dir = pkg_resources.files('TCRAFT.references')
trbv_trac_ohset = pd.read_csv(ref_dir.joinpath('TRBV_TRAC_OHset.csv'), index_col='N')
trbc_trav_ohset = pd.read_csv(ref_dir.joinpath('TRBC_TRAV_OHset.csv'), index_col='N')
trbv_trac_vectors = pd.read_csv(ref_dir.joinpath('TRBV_vector_inserts.csv'))
trbc_trav_vectors = pd.read_csv(ref_dir.joinpath('TRAV_vector_inserts.csv'))
trj_seqs = pd.read_csv(ref_dir.joinpath('TRJ_seqs.csv'), index_col='TRJ').to_dict()['seq']
trav_seqs = parse_and_store_fasta(ref_dir.joinpath('IMGT_TRAV_download.fasta'), 'TRAV')
trbv_seqs = parse_and_store_fasta(ref_dir.joinpath('IMGT_TRBV_download.fasta'), 'TRBV')

with open(ref_dir.joinpath('TRAC_protein.txt'), 'r') as infile:
    trac_protein = infile.read().strip()
with open(ref_dir.joinpath('TRAC_nt_conserved.txt'), 'r') as infile:
    trac_nt_conserved = infile.read().strip()
with open(ref_dir.joinpath('TRBC2_protein.txt'), 'r') as infile:
    trbc2_protein = infile.read().strip()
with open(ref_dir.joinpath('P2A.txt'), 'r') as infile:
    p2a_protein = infile.read().strip()
    p2a_protein = str(Seq(p2a_protein).translate())

# Restriction sites to check for in the assembled TCR sequences. None should be present.
all_re_sites = [BsmBI.site, BsaI.site, BbsI.site, SapI.site, XbaI.site, XhoI.site, EcoRI.site, EcoRV.site, NotI.site]


def rev_comp(seq):
    """
    Returns the reverse complement of a given DNA sequence.
    """
    mapper = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([mapper[nt] for nt in seq[::-1]])


def golden_gate_assembly(vector, insert, rec_site, left_oh, right_oh, space):
    """
    Simulates the Golden Gate assembly process by finding the appropriate
    Type IIs restriction sites and ligating the insert sequence into the vector sequence.

    params
        vector: The vector sequence.
        insert: The insert sequence.
        rec_site: The Type IIS restriction site used for assembly.
        left_oh: The left overhang sequence.
        right_oh: The right overhang sequence.
        space: The number of offset nucleotides between the cut site (overhang) and the restriction site
        for this particular enzyme.

    return
        The assembled sequence after Golden Gate Assembly.
    """
    vector = vector.upper()
    insert = insert.upper()
    rec_site = rec_site.upper()
    left_oh = left_oh.upper()
    right_oh = right_oh.upper()
    spacer = '.'*space
    
    vector_left_idx = re.search(f'{left_oh}{spacer}{rev_comp(rec_site)}', vector).span()[0] #include OH
    insert_left_idx = re.search(f'{rec_site}{spacer}{left_oh}', insert).span()[1]
    vector_right_idx = re.search(f'{rec_site}{spacer}{right_oh}', vector).span()[1]
    insert_right_idx = re.search(f'{right_oh}{spacer}{rev_comp(rec_site)}', insert).span()[0] #include OH
    return vector[:vector_left_idx+4] + insert[insert_left_idx:insert_right_idx+4] + vector[vector_right_idx:]


def prepare_native_trbv_seq(trbv_name): 
    """
    Prepares native TRBV sequence up to and including the 'C' in the CDR3
    params
        trbv_name: ID of the TRBV allele.

    return
        TRBV sequence up to and including the 'C' in the CDR3.
    """
    trbv_aa = trbv_seqs[trbv_name]
    cdr3_start_pos = trbv_aa.rfind('C') + 1
    trimmed_seq = trbv_aa[:cdr3_start_pos]
    return trimmed_seq


def prepare_native_trav_seq(trav_name):
    """
    Prepares native TRAV sequence up to and including the 'C' in the CDR3
    params
        trav_name: ID of the TRBV allele.

    return
        TRAV sequence up to and including the 'C' in the CDR3.
    """
    trav_aa = trav_seqs[trav_name]
    cdr3_start_pos = trav_aa.rfind('C') + 1
    trimmed_seq = trav_aa[:cdr3_start_pos]
    return trimmed_seq


def check_re_sites(seq, sites):
    """
    Checks if any of the given restriction sites are present in the sequence.
    params
        seq: The sequence to check.
        sites: A list of restriction sites to check for.

    return
        A tuple containing the position of the first found restriction site and the site itself.
        If no sites are found, returns (-1, None).
    """
    seq = seq.upper()
    for site in sites:
        site = site.upper()
        if site in seq or rev_comp(site) in seq:
            try:
                pos = seq.index(site)
            except ValueError:
                pos = seq.index(rev_comp(site))
            return pos, site
    return -1, None


def main():
    """
    Main function to validate the assembly of TCR sequences.
    """
    #PARSE INPUT ARGUMENTS
    def csv_file(path):
        if not path.lower().endswith('.csv'):
            raise argparse.ArgumentTypeError("Input CSV must have a .csv extension")
        return path

    parser = argparse.ArgumentParser(description='TCRAFT-validate: Validate that TCRAFT CDR3 oligos are error-free and assemble via Golden Gate into the desired TCR sequence.')
    parser.add_argument('input_csv', type=csv_file, help='CSV file generated from TCRAFT-generate oligo script to validate')
    date_str = pd.Timestamp.now().strftime('%Y%m%d')
    parser.add_argument('--output_dir', type=str, help='Directory to save output data. Default is ./TCRAFT_output_<date>', default=f'./TCRAFT-validate_{date_str}')
    args = parser.parse_args()
    cdr3_oligos = pd.read_csv(args.input_csv)

    #check if required columns are present
    required_columns = ['TRAV', 'TRBV', 'CDR3A', 'CDR3B', 'TRAJ', 'TRBJ', 'Sequence', 'Pool']
    for col in required_columns:
        if col not in cdr3_oligos.columns:
            raise ValueError(f'Input CSV must contain the following columns: {", ".join(required_columns)}')

    print('*******************')
    print('* TCRAFT-validate *')
    print('*******************')
    print()
    print('Beginning assembly validation...')
    correctly_assembled = 0
    step1a_lengths = []
    step1b_lengths = []
    step2a_lengths = []
    step2b_lengths = []
    TCR_vectors = []
    failed_assemblies = []

    for i in tqdm(range(cdr3_oligos.shape[0])):
        params = cdr3_oligos.iloc[i, :]
        cdr3_oligo = params['Sequence']

        #Step 1 - BbsI TRBV TRAC ligation
        trbv_trac_vector = trbv_trac_vectors['Sequence'].loc[trbv_trac_vectors['TRBV'] == params['TRBV']].values[0] + trac_nt_conserved
        left_oh = trbv_trac_ohset['TRBV'].loc[trbv_trac_ohset['TRBV_name'] == params['TRBV']].values[0]
        right_oh = trbv_trac_ohset['TRAC'].loc[trbv_trac_ohset['TRBV_name'] == params['TRBV']].values[0]
        step1_pool = trbv_trac_ohset['set'].loc[trbv_trac_ohset['TRBV_name'] == params['TRBV']].values[0]
        step1_vector = golden_gate_assembly(trbv_trac_vector, cdr3_oligo, BbsI.site, left_oh, right_oh, 2)

        #Step 2 - BsaI/BsmBI TRBC TRAV ligation
        trbc_trav_insert = trbc_trav_vectors['Sequence'].loc[trbc_trav_vectors['TRAV'] == params['TRAV']].values[0]
        left_oh = trbc_trav_ohset['TRBC'].loc[trbc_trav_ohset['TRAV_name'] == params['TRAV']].values[0]
        right_oh = trbc_trav_ohset['TRAV'].loc[trbc_trav_ohset['TRAV_name'] == params['TRAV']].values[0]
        step2_pool = trbc_trav_ohset['set'].loc[trbc_trav_ohset['TRAV_name'] == params['TRAV']].values[0]
        
        if step2_pool == 'A':
            rec_site = BsmBI.site
        else:
            rec_site = BsaI.site

        step2_vector = golden_gate_assembly(step1_vector, trbc_trav_insert, rec_site, left_oh, right_oh, 1)
        gg_tcr_assembly = str(Seq(step2_vector).translate())

        #Validation: Build Reference TCR from building blocks.
        trbv_seq = prepare_native_trbv_seq(params['TRBV'])
        trav_seq = prepare_native_trav_seq(params['TRAV'])
        TCRbeta = trbv_seq + params['CDR3B'][1:] + trj_seqs[params['TRBJ']] + trbc2_protein
        TCRalpha = trav_seq + params['CDR3A'][1:] + trj_seqs[params['TRAJ']] + trac_protein
        real_TCR_assembly = TCRbeta + p2a_protein + TCRalpha + '*' #add beta chain, P2A, alpha chain, and stop codon

        assembly_success = True
        issue = ""
        try:
            assert gg_tcr_assembly == real_TCR_assembly
        except AssertionError:
            print(f'Assembly Failed for TCR # {i+1}:')
            print('Assembled TCR protein sequence:')
            print(gg_tcr_assembly)
            print('Desired TCR protein sequence:')
            print(real_TCR_assembly)
            assembly_success = False
            issue += 'Mismatch with reference TCR sequence;'

        try:
            assert check_re_sites(step2_vector, all_re_sites)[0] == -1
        except AssertionError:
            print(f'Extra restriction site found for TCR # {i+1}:')
            print(step2_vector)
            print('Restriction Site Analysis:')
            print(check_re_sites(step2_vector, all_re_sites))
            assembly_success = False
            issue += 'Extra restriction site found'
        
        if assembly_success: #if no assertion errors, then the TCR was correctly assembled!
            correctly_assembled += 1
            TCR_vector_entry = params.copy()
            TCR_vector_entry.drop('Sequence', inplace=True)
            TCR_vector_entry['Full_TCR_Sequence'] = step2_vector
            TCR_vectors.append(TCR_vector_entry)

            if step1_pool == 'A':
                step1a_lengths.append(len(step1_vector))
            else:
                step1b_lengths.append(len(step1_vector))

            if step2_pool == 'A':
                step2a_lengths.append(len(step2_vector))
            else:
                step2b_lengths.append(len(step2_vector))
            
        else:
            TCR_vector_entry = params.copy()
            TCR_vector_entry.drop('Sequence', inplace=True)
            TCR_vector_entry['Issue'] = issue
            TCR_vector_entry['Oligo_Sequence'] = cdr3_oligo
            TCR_vector_entry['Assembled_TCR_protein'] = gg_tcr_assembly
            TCR_vector_entry['Reference_TCR_protein'] = real_TCR_assembly
            TCR_vector_entry['Assembled_TCR_Sequence'] = step2_vector
            failed_assemblies.append(TCR_vector_entry)

    print('Completed assembly validation.')
    print(f'{correctly_assembled} out of {i+1} TCRs correctly assembled.')
    
    TCR_vectors = pd.DataFrame(TCR_vectors)
    failed_assemblies = pd.DataFrame(failed_assemblies)

    #SAVE OUTPUTS
    print(f'Saving output files to {args.output_dir}')
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        pass
    TCR_vectors.to_csv(os.path.join(args.output_dir, 'Reference_Assembled_TCR_Sequences.csv'), index=False)
    failed_assemblies.to_csv(os.path.join(args.output_dir, 'Failed_Assemblies.csv'), index=False)


    with open(os.path.join(args.output_dir, 'Assembly_Metadata.csv'), 'w') as out:
        out.write('Entry Name,Value\n')
    
        out.write(f'Num Vectors: Step 1A Product,{len(step1a_lengths)}\n')
        out.write(f'Num Vectors: Step 1B Product,{len(step1b_lengths)}\n')
        out.write(f'Avg TRBV-CDR3-TRAC Size: Step 1A Product,{np.mean(step1a_lengths)}\n')
        out.write(f'Avg TRBV-CDR3-TRAC Size: Step 1B Product,{np.mean(step1b_lengths)}\n')

        out.write(f'Num Vectors: Step 2A Product,{len(step2a_lengths)}\n')
        out.write(f'Num Vectors: Step 2B Product,{len(step2b_lengths)}\n')
        out.write(f'Avg TCRb-P2A-TRAC Size: Step 2A Product,{np.mean(step2a_lengths)}\n')
        out.write(f'Avg TCRb-P2A-TRAC Size: Step 2B Product,{np.mean(step2b_lengths)}\n')

if __name__ == '__main__':
    main()
