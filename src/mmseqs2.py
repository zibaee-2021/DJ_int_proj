import os, subprocess, glob
import shutil
from typing import List, Tuple
from collections import defaultdict
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from parso.python.tree import WithStmt

import pdb_model_stats as pms

# NOTE I HAVE ANOTHER DICT FOR SOME HETATM WHICH ARE POST-TRANSLATIONALLY-MODIFIED RESIDUES.
# E.G. I CONVERT 'PTR' (PHOSPHOTYROSINE) TO 'Y' RATHER THAN LEAVING IT OUT:
# "Pyroglutamic acid is almost always derived from N-terminal glutamine."
# ON THE OTHER HAND, THERE ARE RESIDUES THAT ARE MAN-MADE SYNTHETIC RESIDUES.
# I HAVE A DICT FOR THESE AS WELL MAPPING TO NATURAL RESIDUES, BUT FOR NOW I AM OPTING TO EXCLUDE THE WHOLE CHAIN.
def aa3to1_unused(pdbid: str, chain: str, three_char_seq: List[Tuple[int, str]]):
    aa_3to1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
        'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
        'TRP': 'W', 'TYR': 'Y'}

    aa_3to1_natural_ptm = {
        'ASQ': 'S',  # Acetylserine
        'CGU': 'E',  # Gamma-carboxyglutamate
        'CIR': 'R',  # Citrulline
        'CSO': 'C',  # S-hydroxycysteine
        'CSX': 'C',  # Cysteine S-oxide
        'FME': 'M',  # N-formylmethionine
        'HIC': 'H',  # 4-methylhistidine
        'HYP': 'P',  # Hydroxyproline
        'KCX': 'K',  # N-epsilon-carboxylysine
        'MLY': 'K',  # Alternate for methyllysine
        'MLZ': 'K',  # N6-methyllysine
        'MSE': 'M',  # Selenomethionine
        'PCA': 'Q',  # Pyroglutamic acid
        'PTR': 'Y',  # Phosphotyrosine
        'PYR': 'Q',  # Alternate code for pyroglutamic acid
        'SEP': 'S',  # Phosphoserine
        'TPO': 'T',  # Phosphothreonine
    }

    aa_3to1_manmade = {
        'DAL': 'A',  # D-alanine
        'DAR': 'R',  # D-arginine
        'DCY': 'C',  # D-cysteine
        'DGL': 'E',  # D-glutamic acid
        'DGN': 'Q',  # D-glutamine
        'DHI': 'H',  # D-histidine
        'DIL': 'I',  # D-isoleucine
        'DLE': 'L',  # D-leucine
        'DLY': 'K',  # D-lysine
        'DPH': 'F',  # D-phenylalanine
        'DSN': 'S',  # D-serine
        'DTH': 'T',  # D-threonine
        'DTY': 'Y',  # D-tyrosine
        'DVA': 'V',  # D-valine
    }

    fasta_seq = []
    for i, (_, res) in enumerate(three_char_seq):
        if res not in aa_3to1:
            # raise KeyError(f'{res} not found in list of three-letter codes. PDBid={pdbid}, Chain={cif_chain}')
            print(f'{pdbid}_{chain}: {res} not one of the 20 natural residues.')
            if res not in aa_3to1_natural_ptm:
                print(f'{pdbid}_{chain}: {res} also not in my list of 17 post-translationally modified residues.')
                if res not in aa_3to1_manmade: # not natural and not PTM - so remove:
                    print(f'{pdbid}_{chain}: {res} also not in my list of 14 synthetic residues. '
                          f'Therefore this pdbid_chain will be removed, as I have no idea what it is !')
                    return None
                else:
                    print(f'{pdbid}_{chain}: {res} is in my list of 14 synthetic residues.')
                    print(f'Therefore this pdbid_chain will be removed, as it contains synthetic residues.')
                    return None
            else: # natural PTM residue:
                aa1char = aa_3to1_natural_ptm[res]
                print(f'{pdbid}_{chain}: {res} is in PTM residues list and is being mapped to {aa1char}')
                if (i + 1) % 80 == 0:
                    fasta_seq.append(aa1char + '\n')
                else:
                    fasta_seq.append(aa1char)
        else:  # natural residue:
            if (i + 1) % 80 == 0:
                fasta_seq.append(aa_3to1[res] + '\n')
            else:
                fasta_seq.append(aa_3to1[res])

    return ''.join(fasta_seq)


def map_aa3to1(aa3seq: list, pdbid: str, chain: str) -> str:
    aa3to1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
        'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
        'TRP': 'W', 'TYR': 'Y'}
    fasta_seq = []
    for aa3 in aa3seq:
        if aa3 not in aa3to1:
            raise KeyError(f'{aa3} not found in list of three-letter codes. PDBid={pdbid}, Chain={chain}. '
                           f'So is it not one of the 20 natural residues ??')
        else:
            if fasta_seq and len(fasta_seq) % 80 == 0:
                fasta_seq.append(aa3to1[aa3] + '\n')
            else:
                fasta_seq.append(aa3to1[aa3])
    return ''.join(fasta_seq)



def rp_rawcifs_dir(het_hom: str) -> str:
    return os.path.join('..', 'data', 'NMR', 'raw_cifs', het_hom)


def rp_tokcifs_dir(het_hom: str) -> str:
    return os.path.join('..', 'data', 'NMR', 'tokenised_cifs', het_hom)


def rp_mmseqs_dir(het_hom: str) -> str:
    return os.path.join('..', 'data', 'NMR', 'mmseqs', het_hom)


def rp_mmseqs_fasta_dir(het_hom: str) -> str:
    return os.path.join(rp_mmseqs_dir(het_hom), 'fasta')


def write_fasta(het_hom: str, pdbid: str, chains: list) -> str:
    toks_dir = rp_tokcifs_dir(het_hom)
    fasta_seq = []
    for chain in chains:
        pdbid_chain_pdf = pd.read_csv(os.path.join(toks_dir, f'{pdbid}_{chain}.ssv'), sep=' ')
        models = pdbid_chain_pdf['A_pdbx_PDB_model_num'].unique().tolist()
        pdf_one_model = pdbid_chain_pdf[pdbid_chain_pdf['A_pdbx_PDB_model_num'] == models[0]]
        aa3seq = pdf_one_model['S_mon_id'].tolist()
        aa1seq = map_aa3to1(aa3seq, pdbid, chain)
        fasta_seq.append(f'>{pdbid}_{chain}')
        fasta_seq.append(aa1seq)

    rp_fasta_file = os.path.join(rp_mmseqs_fasta_dir(het_hom), f'{pdbid}.fasta')

    with open(rp_fasta_file, 'w') as f:
        f.writelines('\n'.join(fasta_seq) + '\n')

    return rp_fasta_file


def run_mmseqs_all_vs_all(rp_fasta_f, pdbid: str):
    """
    Runs MMseqs2 all-vs-all search on the provided FASTA file.
    Outputs:
      - results.m8 (tabular results)
      - Pandas DataFrame of alignments
    """
    het_hom = rp_fasta_f.split('/')[4]
    rp_dir2delete = os.path.join(rp_mmseqs_dir(het_hom), 'dir2delete')
    shutil.rmtree(rp_dir2delete, ignore_errors=True)

    os.makedirs(rp_dir2delete, exist_ok=True)
    tmp_dir = os.path.join(rp_dir2delete, f'{pdbid}_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    db_name = os.path.join(rp_dir2delete, f'{pdbid}_db')
    result_db = os.path.join(rp_dir2delete, f'{pdbid}_result')
    output_m8 = os.path.join(rp_dir2delete, f'{pdbid}_results.m8')

    # 1. Create database:
    subprocess.run(['mmseqs', 'createdb',
                    rp_fasta_f, db_name],
                   stdout=subprocess.DEVNULL, check=True)

    # 2. Search: all-vs-all:
    # subprocess.run(['mmseqs', 'search',
    #     db_name,   # query db
    #     db_name,   # target db
    #     result_db, tmp_dir, '--threads', '4'],
    #                stdout=subprocess.DEVNULL, check=True)

    # 2. Search: all-vs-all (with lower threshold):
    subprocess.run(['mmseqs', 'search',
                    db_name, db_name,
                    result_db, tmp_dir,
                    '--alignment-mode', '3',
                    '--threads', '4',
                    '-s', '7.5',
                    '-e', '1000'],
                   stdout=subprocess.DEVNULL, check=True)

    # 3. Convert alignments to tabular output
    subprocess.run(['mmseqs', 'convertalis',
                    db_name, db_name,
                    result_db, output_m8,
                    '--format-output',
                    'query,target,evalue,pident,alnlen'],
                   stdout=subprocess.DEVNULL, check=True)

    # 4. Read results into Pandas
    df = pd.read_csv(output_m8, sep='\t', header=None, names=['query','target','evalue','pident','alnlen'])
    return df


def filter_and_write_results(pdbid: str, het_hom: str, pdf):
    pdf = pdf[pdf['query'] != pdf['target']]
    rp_mmseqs_dir_ = rp_mmseqs_dir(het_hom)
    if not pdf.empty:
        results_dir = os.path.join(rp_mmseqs_dir_, 'results')
        os.makedirs(results_dir, exist_ok=True)
        pdf.to_csv(os.path.join(results_dir, f'{pdbid}.csv'), index=False)
    else:
        print(f'No results for {pdbid}. Adding id to list file.')
        rp_zero_idty_lst = os.path.join(rp_mmseqs_dir_, 'PDBid_no_idty.lst')

        if os.path.exists(rp_zero_idty_lst):
            with open(rp_zero_idty_lst, 'r') as f:
                pdbids_no_idty = f.readlines()
            pdbids_no_idty = [_pdbid.removesuffix('\n') for _pdbid in pdbids_no_idty]
            pdbids_no_idty.append(pdbid)
            pdbids_no_idty.sort()
            with open(rp_zero_idty_lst, 'w') as f:
                f.write('\n'.join(pdbids_no_idty) + '\n')
        else:
            with open(rp_zero_idty_lst, 'w') as f:
                f.write(f'{pdbid}\n')
    return pdf


if __name__ == '__main__':
    write_fasta(het_hom='heteromeric', pdbid='1A0N', chains=['A', 'B'])
    # het_dir = os.path.join('..', 'data', 'NMR', 'mmseqs', 'heteromeric')
    # lst_f = os.path.join(het_dir, 'PDBid_no_idty.lst')
    # if os.path.exists(lst_f):
    #     os.remove(lst_f)
    #     print(f'{os.path.basename(lst_f)} deleted')
    # else:
    #     print(f'{os.path.basename(lst_f)} not found')
    #
    # fasta_dir = os.path.join(het_dir, 'fasta')
    # if os.path.isdir(fasta_dir):
    #     shutil.rmtree(fasta_dir)
    #     print(f'{os.path.basename(fasta_dir)} deleted.')
    # else:
    #     print(f'{os.path.basename(fasta_dir)} does not exist.')
    #
    # rp_pdbid_chains = os.path.join('..', 'data', 'NMR', 'multimodel_lists', 'het_multimod_2104_pdbid_chains.txt')
    # het_hom = os.path.basename(rp_pdbid_chains)[:3]
    # het_hom = 'heteromeric' if het_hom == 'het' else 'homomeric'
    #
    # pdbid_chains_dict = pms.pdbid_dict_chain(rp_pdbid_chains)
    # # 2104 heteromeric pdbid_chains --> 932 heteromeric pdbids.
    # pdbids = list(pdbid_chains_dict.keys())
    #
    # print(f'Starting to process {len(pdbids)} PDBids to write to fasta.')
    # for i, _pdbid in enumerate(pdbids):
    #     print(f'Processing {i}: {_pdbid}')
    #     # MMseqs2
    #     rp_fasta_file = write_fasta(het_hom, _pdbid, pdbid_chains_dict)
    #     # result_pdf = run_mmseqs_all_vs_all(rp_fasta_file, _pdbid)
    #     # filtered_pdf = filter_and_write_results(_pdbid, het_hom, result_pdf)
    #
    #     # for chain in pdbid_chains_dict[_pdbid]:
    #     #     if not filtered_pdf.empty:
    #     #         identity = 'see idnty csv'
    #     #     else:
    #     #         identity = '0'



#     _het_hom = 'heteromeric'
#     # _pdbid = '1A0N'
#     # _pdbid_chains_dict ={_pdbid: ['A', 'B']}
#     _pdbid = '1AOU'
#     _pdbid_chains_dict = {_pdbid: ['A', 'B']}
#     _rp_fasta_file = write_fasta(_het_hom, _pdbid, _pdbid_chains_dict)
#     pass
    # df_results = run_mmseqs_all_vs_all(rp_fasta_f=__rp_fasta_file, pdbid=_pdbid)
    # print(df_results)
    # filter_results(_pdbid, _het_hom, df_results)
    # df_results = df_results[df_results['query'] != df_results['target']]
    #
    # _rp_mmseqs_dir = relpath_mmseqs_dir(_het_hom)
    # if not df_results.empty:
    #     _pdf_csv_results_dir = os.path.join(_rp_mmseqs_dir, 'results')
    #     os.makedirs(_pdf_csv_results_dir, exist_ok=True)
    #     df_results.to_csv(os.path.join(_pdf_csv_results_dir, f'{_pdbid}.csv'), index=False)
    # else:
    #     print(f'No results for {_pdbid}. Adding id to list file.')
    #     _rp_zero_idty_lst = os.path.join(_rp_mmseqs_dir, 'PDBid_no_idty.lst')
    #
    #     if os.path.exists(_rp_zero_idty_lst):
    #         with open(_rp_zero_idty_lst, 'r') as f:
    #             _pdbids_no_idty = f.readlines()
    #         _pdbids_no_idty = [__pdbid.removesuffix('\n') for __pdbid in _pdbids_no_idty]
    #         _pdbids_no_idty.append(_pdbid)
    #         _pdbids_no_idty.sort()
    #         with open(_rp_zero_idty_lst, 'w') as f:
    #             f.write('\n'.join(_pdbids_no_idty) + '\n')
    #     else:
    #         with open(_rp_zero_idty_lst, 'w') as f:
    #             f.write(f'{_pdbid}\n')


    # NOTE: Using Biopython results in use of a chain name attribute of _atom_site which seems problematic.
    # Specifically _auth_asym_id instead of _label_asym_id (the latter is what I used in my project).
    # So the following 4 lines were removed from write_fasta(). Instead I use MMCIF2DICT.
    # parser = MMCIFParser(QUIET=True)
    # bio_struct = parser.get_structure('', rp_cif)
    # model = bio_struct[0]
    # for chain in model: ...
