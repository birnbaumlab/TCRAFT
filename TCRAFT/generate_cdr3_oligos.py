#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#IMPORTS AND FILE LOADING
import numpy as np
import random
import pandas as pd
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.Restriction import *
import json
from Bio.SeqUtils import gc_fraction
from tqdm import tqdm
import argparse
import os
import importlib.resources as pkg_resources

"""
This section loads the reference files needed for oligo assembly.
"""
ref_dir = pkg_resources.files('TCRAFT.references')
with open(ref_dir.joinpath('H_sapiens_codon_usage.json'), 'r') as infile:
    codon_freq_dict = json.load(infile)

trbv_trac_ohset = pd.read_csv(ref_dir.joinpath('TRBV_TRAC_OHset.csv'), index_col='N')
trbc_trav_ohset = pd.read_csv(ref_dir.joinpath('TRBC_TRAV_OHset.csv'), index_col='N')
trbv_trac_OH_frags = pd.read_csv(ref_dir.joinpath('TRBV_TRAC_OH_frags.csv'))
trbc_trav_OH_frags = pd.read_csv(ref_dir.joinpath('TRBC_TRAV_OH_frags.csv'))
trj_seqs = pd.read_csv(ref_dir.joinpath('TRJ_seqs.csv'))

with open(ref_dir.joinpath('Oligo_Components.json'), 'r') as infile:
    oligo_components = json.load(infile)

#Load primer and restriction site sequences
fwd_primer = oligo_components['fwd_primer']
rev_primer = oligo_components['rev_primer']
bbs1 = oligo_components['bbs1']
bsa1 = oligo_components['bsa1']
bsmb1 = oligo_components['bsmb1']


#HELPER FUNCTIONS
def rev_comp(seq):
    """
    Returns the reverse complement of a DNA sequence 'seq'.
    """
    mapper = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([mapper[nt] for nt in seq[::-1]])

def dfs_seq_builder(ref_aa, primer='', aa_pos=0):
    """
    This function uses a depth-first search to build all possible DNA sequences that translate to the given amino acid sequence.

    params:
        ref_aa: str, the reference amino acid sequence
        primer: str, the current DNA sequence being built
        aa_pos: int, the current position in the amino acid sequence

    returns:
        list of all possible DNA sequences that translate to the given amino acid sequence
    """
    if aa_pos >= len(ref_aa):
        return [primer]
    output = []
    for codon_option, _ in codon_freq_dict[ref_aa[aa_pos]]:
        output.extend(dfs_seq_builder(ref_aa, primer+codon_option, aa_pos+1))
    return output 


#FUNCTIONS FOR CDR3 CONSTRUCTION
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
    for site in sites:
        if site in seq or rev_comp(site) in seq:
            try:
                pos = seq.index(site)
            except ValueError:
                pos = seq.index(rev_comp(site))
            return pos, site
    return -1, None


def check_homopolymer_sites(seq):
    """
    Checks if homopolymers >= 5bp in length are present in the sequence.
    params
        seq: The sequence to check.

    return
        A tuple containing the position of the first found homopolymer site and the site itself.
        If no sites are found, returns (-1, None).
    """
    for site in ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']:
        if site in seq:
            try:
                pos = seq.index(site)
            except ValueError:
                pos = seq.index(rev_comp(site))
            return pos, site
    return -1, None


def remove_restriction_sites(seq, aa, sites):
    """
    Iteratively removes restriction sites from the sequence by replacing them with a random codon that does not contain the restriction site.
    params:
        seq: The sequence to check.
        aa: The amino acid sequence corresponding to the DNA sequence.
        sites: A list of restriction sites to check for and remove.

    return:
        The modified sequence with restriction sites removed.
    """
    iters = 0
    while check_re_sites(seq, sites)[0] != -1:
        pos, site = check_re_sites(seq, sites)
        start = pos // 3
        end = (pos + len(site)) // 3
        
        replace = None
        for option in dfs_seq_builder(aa[start:end+1]):
            start_minus_one = max(0, start-1)
            end_plus_one = min(end+1, len(seq)//3 - 1)
            option_with_context = seq[start_minus_one*3:start*3] + option + seq[end*3:end_plus_one*3]
            if check_re_sites(option_with_context, sites)[0] == -1:
                replace = option

        if iters >= 100 or replace is None:
            raise ValueError(f'Note: Failed to remove restriction site. Seq = {seq}.')
        
        seq = seq[:start*3] + replace + seq[(end+1)*3:]
        iters += 1
        
    assert Seq(seq).translate() == aa
    return seq


def remove_homopolymer_sites(seq, aa, sites):
    """
    Iteratively removes homopolymer sites from the sequence by replacing them with a random codon that does not contain the homopolymer site.
    Simultaneously checks for restriction sites to ensure they are not reintroduced after removing homopolymers.

    params:
        seq: The sequence to check.
        aa: The amino acid sequence corresponding to the DNA sequence.
        sites: A list of restriction sites to check for and remove. This is to ensure that the homopolymer replacement does not introduce new restriction sites.

    return:
        The modified sequence with restriction sites removed.
    """
    iters = 0
    while check_homopolymer_sites(seq)[0] != -1:
        pos, site = check_homopolymer_sites(seq)
        start = pos // 3
        end = (pos + len(site)) // 3
        
        replace = None
        for option in dfs_seq_builder(aa[start:end+1]):
            start_minus_one = max(0, start-1)
            end_plus_one = min(end+1, len(seq)//3 - 1)
            option_with_context = seq[start_minus_one*3:start*3] + option + seq[end*3:end_plus_one*3]
            if check_re_sites(option_with_context, sites)[0] == -1 and check_homopolymer_sites(option_with_context)[0] == -1:
                replace = option

        if iters >= 100 or replace is None:
            raise ValueError(f'Note: Failed to remove homopolymer. Seq = {seq}')
        
        seq = seq[:start*3] + replace + seq[(end+1)*3:]
        iters += 1
        
    assert Seq(seq).translate() == aa
    return seq


def generate_cleaned_DNA(aa_seq, freq_cutoff=0.15):
    """
    Generates a cleaned DNA sequence from the given amino acid sequence.
    The DNA sequence is generated by sampling from the codon usage table, ensuring that the GC content is between 25% and 65%.
    The generated sequence is then cleaned by removing restriction sites and homopolymer sites.
    params:
        aa_seq: str, the amino acid sequence to generate DNA for.
        freq_cutoff: float, the minimum frequency of a codon to be considered for sampling.

    return:
        The cleaned DNA sequence as a string
    """
    
    def generate_DNA(): #helper function to generate a DNA sequence for the given amino acid sequence.
        seq = ""
        while seq == "" or gc_fraction(seq) < 0.25 or gc_fraction(seq) > 0.65:
            seq = ""
            for aa in aa_seq:
                codons = []
                freqs = []
                for codon, freq in codon_freq_dict[aa]:
                    if freq > freq_cutoff:
                        codons.append(codon)
                        freqs.append(freq)
                freqs = (np.array(freqs) / np.sum(freqs)).tolist()
                seq += random.choices(population=codons, weights=freqs, k=1)[0]
        assert Seq(seq).translate() == aa_seq
        return seq

    #clean the DNA
    tries = 0
    while True:
        try:
            seq = generate_DNA()
            sites = [BsaI.site, BsmBI.site, BbsI.site, SapI.site, XbaI.site, EcoRI.site, XhoI.site, EcoRV.site, NotI.site]
            seq_minus_re = remove_restriction_sites(seq, aa_seq, sites)
            final_seq = remove_homopolymer_sites(seq_minus_re, aa_seq, sites)
            assert check_re_sites(seq, sites)[0] == -1 and check_homopolymer_sites(seq)[0] == -1
            break
        except Exception as e:
            tries += 1
            if tries > 50:
                print(e)
                raise ValueError('Note: Failed to generate cleaned DNA after 50 sampling attempts.')
    
    return final_seq


def concatenate_oligo(oligo):
    """
    Concatenates the oligo fragments into a single sequence.
    params:
        oligo: dict, the oligo fragments to concatenate.
    return:
        The concatenated oligo sequence as a string.
    """
    final_oligo_seq = fwd_primer + \
                      oligo['left_re'] + \
                      oligo['trbv_right_frag'] + \
                      oligo['cdr3bj_seq'] + \
                      oligo['trbc_left_frag'] + \
                      oligo['middle_re'] + \
                      oligo['trav_right_frag'] + \
                      oligo['cdr3aj_seq'] + \
                      oligo['trac_left_frag'] + \
                      oligo['right_re'] + \
                      rev_primer
    return final_oligo_seq


def test_oligo(oligo_seq, trav_set):
    """
    Tests the fully assembled oligo sequence for restriction sites.
    params:
        oligo_seq: str, the oligo sequence to test.
        trav_set: str, the TRAV set to check for the right restriction sites.
    return:
        True if the oligo sequence passes the tests, False otherwise.
    """
    # check restriction sites
    rb = BsmBI + BsaI + BbsI + SapI + XbaI + EcoRI + XhoI + NotI + EcoRV
    restriction_analysis = Analysis(rb, Seq(oligo_seq)).full()
    
    try:
        for key in restriction_analysis.keys():
            if key == BbsI: #need to patch a structural problem with trav12-3
                assert len(restriction_analysis[key]) == 2
            elif key == BsaI:
                if trav_set == 'A':
                    assert len(restriction_analysis[key]) == 1
                else: 
                    assert len(restriction_analysis[key]) == 2
            elif key == BsmBI:
                if trav_set == 'A':
                    assert len(restriction_analysis[key]) == 2
                else: 
                    assert len(restriction_analysis[key]) == 1
            else: #all others minus XhoI, XbaI
                assert len(restriction_analysis[key]) == 0
    except AssertionError:
        return False
    return True


def create_CDR3_oligo(trav_name, trbv_name, cdr3a, cdr3b, traj_name, trbj_name):
    """
    Creates a CDR3 oligo from the given TCR information.
    params:
        trav_name: TRAV allele ID
        trbv_name: TRBV allele ID
        cdr3a: CDR3 alpha amino acid sequence, starting from C and ending with F
        cdr3b: CDR3 beta amino acid sequence, starting from C and ending with F
        traj_name: TRAJ allele ID
        trbj_name: TRBJ allele ID
    return:
        A tuple containing fully assembled oligo sequence as a string and the TRBV set for Step 1 pooling
    """
    oligo_frags = {}

    if trav_name not in trbc_trav_OH_frags['trav'].values or trbv_name not in trbv_trac_OH_frags['trbv'].values:
        raise ValueError(f'Note: Invalid TRAV {trav_name} and/or TRBV {trbv_name}, skipping this TCR.')
    
    oligo_frags['trbv_right_frag'] = trbv_trac_OH_frags['trbv_right_frag'].loc[trbv_trac_OH_frags['trbv'] == trbv_name].values[0]
    oligo_frags['trbc_left_frag'] = trbc_trav_OH_frags['trbc_left_frag'].loc[trbc_trav_OH_frags['trav'] == trav_name].values[0]
    oligo_frags['trav_right_frag'] = trbc_trav_OH_frags['trav_right_frag'].loc[trbc_trav_OH_frags['trav'] == trav_name].values[0]
    oligo_frags['trac_left_frag'] = trbv_trac_OH_frags['trac_left_frag'].loc[trbv_trac_OH_frags['trbv'] == trbv_name].values[0]

    #Find TRAJ and TRBJ
    try:
        traj = trj_seqs.loc[trj_seqs.TRJ == traj_name, 'seq'].values[0]
        trbj = trj_seqs.loc[trj_seqs.TRJ == trbj_name, 'seq'].values[0]

    except IndexError:
        raise ValueError(f'Note: Invalid TRAJ {traj_name} and/or TRBJ {trbj_name}, skipping this TCR.')
    

    #prepare full CDR3A+J and CDR3B+J seq, cut off the first C (already in the TRAV/TRBV)
    #residue between TRAJ and TRAC is specified in the TRAJ
    oligo_frags['cdr3aj_seq'] = generate_cleaned_DNA(cdr3a[1:] + traj).lower()
    oligo_frags['cdr3bj_seq'] = generate_cleaned_DNA(cdr3b[1:] + trbj).lower()

    #patch frequent edge cases that occur in the J-C junction to prevent excessive resampling
    if oligo_frags['cdr3bj_seq'][-3:] == 'ctc': #patch frequent CTCGAG edge case in TRBJ-TRBC junction
        oligo_frags['cdr3bj_seq'] = oligo_frags['cdr3bj_seq'][:-3] + 'ctg' #another leucine codon 
    if oligo_frags['cdr3aj_seq'][-3:] == 'gat': #patch frequent GATATC edge case in TRAJ-TRAC junction
        oligo_frags['cdr3aj_seq'] = oligo_frags['cdr3aj_seq'][:-3] + 'gac' #another aspartate codon

    #generate restriction site sequences
    trbv_set = trbv_trac_ohset.loc[trbv_trac_ohset.TRBV_name == trbv_name, 'set'].values[0]
    trav_set = trbc_trav_ohset.loc[trbc_trav_ohset.TRAV_name == trav_name, 'set'].values[0]

    # all CDRs have BbsI for step 1
    oligo_frags['left_re'], oligo_frags['right_re'] = bbs1

    # assign middle restriction site based on TRAV set
    if trav_set == 'A':
        oligo_frags['middle_re'] = bsmb1
    else:
        oligo_frags['middle_re'] = bsa1

    #resample CDR3s if restriction sites are present in the junctions between sequence fragments
    resample_count = 0
    while not test_oligo(concatenate_oligo(oligo_frags), trav_set):
        oligo_frags['cdr3aj_seq'] = generate_cleaned_DNA(cdr3a[1:] + traj).lower()
        oligo_frags['cdr3bj_seq'] = generate_cleaned_DNA(cdr3b[1:] + trbj).lower()
        resample_count += 1
        if resample_count > 50:
            raise ValueError('Note: Failed to fix restriction sites in oligo junctions after 50 resampling attempts. Skipping this TCR.')
        
    if resample_count > 0:
        print(f'Resampled {resample_count} times to fix restriction site issue.')

    return concatenate_oligo(oligo_frags), trbv_set #to split Step1 pools


def main():
    """
    Main function to generate CDR3 oligos from TCR information.
    """
    #PARSE INPUT ARGUMENTS
    def csv_file(path):
        if not path.lower().endswith('.csv'):
            raise argparse.ArgumentTypeError("Input CSV must have a .csv extension")
        return path

    parser = argparse.ArgumentParser(description='TCRAFT-generate: Generate CDR3 Oligos from TCR sequence annotation data.')
    parser.add_argument('input_csv', type=csv_file, help='Input CSV file with TCR information. Must contain the following columns: V_alpha, V_beta, CDR3_alpha, CDR3_beta, J_alpha, J_beta')
    date_str = pd.Timestamp.now().strftime('%Y%m%d')
    parser.add_argument('--output_dir', type=str, help='Directory to save output files. Default is ./TCRAFT_output_<date>', default=f'./TCRAFT-generate_{date_str}')
    args = parser.parse_args()
    tcr_list = pd.read_csv(args.input_csv)

    #check if V_alpha, V_beta, CDR3_alpha, CDR3_beta, J_alpha, J_beta columns are present in tcr_list
    required_columns = ['V_alpha', 'V_beta', 'CDR3_alpha', 'CDR3_beta', 'J_alpha', 'J_beta']
    for col in required_columns:
        if col not in tcr_list.columns:
            raise ValueError(f'Input CSV must contain the following columns: {", ".join(required_columns)}')

    #subset tcr_list to only include required columns
    tcr_list = tcr_list[required_columns]

    #drop rows with NA values
    tcr_list = tcr_list.dropna()

    #MAIN LOOP - GENERATE CDR3 OLIGOS
    print('*******************')
    print('* TCRAFT-generate *')
    print('*******************')
    print()
    print('Generating CDR3 oligos...')
    
    invalid_cdr3_num = 0
    cdr3_oligo_list = {'TRAV':[], 'TRBV':[], 'CDR3A':[], 'CDR3B':[], 'TRAJ':[], 'TRBJ':[], 'Length':[], 'GC':[], 'Pool':[], 'Sequence':[]}
    invalid_cdr3_oligo_list = {'TRAV':[], 'TRBV':[], 'CDR3A':[], 'CDR3B':[], 'TRAJ':[], 'TRBJ':[], 'Issue':[]}

    for i in tqdm(range(tcr_list.shape[0])):
        tcr = tcr_list.iloc[i, :]
        try:
            cdr3_oligo, pool = create_CDR3_oligo(tcr['V_alpha'], tcr['V_beta'], tcr['CDR3_alpha'], tcr['CDR3_beta'], tcr['J_alpha'], tcr['J_beta'])
        
            cdr3_oligo_list['TRAV'].append(tcr['V_alpha'])
            cdr3_oligo_list['TRBV'].append(tcr['V_beta'])
            cdr3_oligo_list['CDR3A'].append(tcr['CDR3_alpha'])
            cdr3_oligo_list['CDR3B'].append(tcr['CDR3_beta'])
            cdr3_oligo_list['TRAJ'].append(tcr['J_alpha'])
            cdr3_oligo_list['TRBJ'].append(tcr['J_beta'])
            cdr3_oligo_list['Length'].append(len(cdr3_oligo))
            cdr3_oligo_list['GC'].append(gc_fraction(cdr3_oligo))
            cdr3_oligo_list['Pool'].append(pool)
            cdr3_oligo_list['Sequence'].append(cdr3_oligo)

        except ValueError as e:
            print(e)
            invalid_cdr3_num += 1

            invalid_cdr3_oligo_list['TRAV'].append(tcr['V_alpha'])
            invalid_cdr3_oligo_list['TRBV'].append(tcr['V_beta'])
            invalid_cdr3_oligo_list['CDR3A'].append(tcr['CDR3_alpha'])
            invalid_cdr3_oligo_list['CDR3B'].append(tcr['CDR3_beta'])
            invalid_cdr3_oligo_list['TRAJ'].append(tcr['J_alpha'])
            invalid_cdr3_oligo_list['TRBJ'].append(tcr['J_beta'])
            invalid_cdr3_oligo_list['Issue'].append(str(e))

    print(f'Generated {len(cdr3_oligo_list["Sequence"])} valid CDR3 oligos.')
    print(f'Failed to generate {invalid_cdr3_num} CDR3 oligos, see log for details.')

    #POST-PROCESSING
    cdr3_oligo_list = pd.DataFrame(cdr3_oligo_list)
    invalid_cdr3_oligo_list = pd.DataFrame(invalid_cdr3_oligo_list)
    cdr3_oligo_list.sort_values('Pool', inplace=True)
    over_300_num = np.sum(pd.DataFrame(cdr3_oligo_list)['Length'] > 300)
    print(f'{over_300_num}/{cdr3_oligo_list.shape[0]} oligos ({round(100*over_300_num/cdr3_oligo_list.shape[0], 2)}%) are over 300 bp long.')

    poolA_under300 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'A', cdr3_oligo_list['Length'] <= 300)]
    poolA_over300 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'A', cdr3_oligo_list['Length'] > 300)]
    poolB_under300 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'B', cdr3_oligo_list['Length'] <= 300)]
    poolB_over300 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'B', cdr3_oligo_list['Length'] > 300)]
    all_under300 = cdr3_oligo_list.loc[cdr3_oligo_list['Length'] <= 300]
    all_over300 = cdr3_oligo_list.loc[cdr3_oligo_list['Length'] > 300]


    #SAVE OUTPUTS
    print(f'Saving output files to {args.output_dir}')
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        pass
    
    cdr3_oligo_list.to_csv(os.path.join(args.output_dir, f'All_{cdr3_oligo_list.shape[0]}_CDR3_oligos.csv'), index=False)
    invalid_cdr3_oligo_list.to_csv(os.path.join(args.output_dir, 'Invalid_CDR3_table.csv'), index=False)
    
    all_under300.to_csv(os.path.join(args.output_dir, f'Under300_All_CDR3_oligos.csv'), index=False)
    all_over300.to_csv(os.path.join(args.output_dir, f'Over300_All_CDR3_oligos.csv'), index=False)

    poolA_under300.to_csv(os.path.join(args.output_dir, 'Under300_PoolA_CDR3_oligos.csv'), index=False)
    poolB_under300.to_csv(os.path.join(args.output_dir, 'Under300_PoolB_CDR3_oligos.csv'), index=False)
    poolA_over300.to_csv(os.path.join(args.output_dir, 'Over300_PoolA_CDR3_oligos.csv'), index=False)
    poolB_over300.to_csv(os.path.join(args.output_dir, 'Over300_PoolB_CDR3_oligos.csv'), index=False)

    #save oligo pool metadata
    with open(os.path.join(args.output_dir, 'Oligo_Pool_Metadata.csv'), 'w') as out:
        out.write('Entry Name,Value\n')

        out.write(f'Num Oligos: Pool A Under 300bp,{poolA_under300.shape[0]}\n')
        out.write(f'Num Oligos: Pool B Under 300bp,{poolB_under300.shape[0]}\n')
        out.write(f'Num Oligos: Pool A Over 300bp,{poolA_over300.shape[0]}\n')
        out.write(f'Num Oligos: Pool B Over 300bp,{poolB_over300.shape[0]}\n')

        out.write(f'Avg Oligo Size: Pool A Under 300bp,{poolA_under300["Length"].mean()}\n')
        out.write(f'Avg Oligo Size: Pool B Under 300bp,{poolB_under300["Length"].mean()}\n')
        out.write(f'Avg Oligo Size: Pool A Over 300bp,{poolA_over300["Length"].mean()}\n')        
        out.write(f'Avg Oligo Size: Pool B Over 300bp,{poolB_over300["Length"].mean()}\n')


if __name__ == "__main__":
    main()
    