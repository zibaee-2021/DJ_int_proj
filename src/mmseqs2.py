"""
Installed mmseqs2 by `conda install -c conda-forge -c bioconda mmseqs2` in activated conda env.
"""
import os, subprocess, glob
from time import time
import shutil
from typing import List, Tuple
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from jupyter_core.migrate import get_ipython_dir
from parso.python.tree import WithStmt
import pdb_model_stats as pms


# NOTE I HAVE ANOTHER DICT FOR SOME HETATM WHICH ARE POST-TRANSLATIONALLY-MODIFIED RESIDUES.
# E.G. I CONVERT 'PTR' (PHOSPHOTYROSINE) TO 'Y' RATHER THAN LEAVING IT OUT:
# "Pyroglutamic acid is almost always derived from N-terminal glutamine."
# ON THE OTHER HAND, THERE ARE RESIDUES THAT ARE MAN-MADE SYNTHETIC RESIDUES.
# I HAVE A DICT FOR THESE AS WELL MAPPING TO NATURAL RESIDUES, BUT FOR NOW I AM OPTING TO EXCLUDE THE WHOLE CHAIN.


def map_aa3to1(aa3seq: list, pdbid: str, chain: str) -> str:
    aa3to1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
        'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
        'TRP': 'W', 'TYR': 'Y'
    }
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

# BUILDING RELATIVE PATHS:
def _rp_nmr_dir():
    return os.path.join('..', 'data', 'NMR')


def _rp_rawcifs_dir(sub_dir) -> str:
    return os.path.join(_rp_nmr_dir(), 'raw_cifs', sub_dir)

def rp_raw_pdbs_dir(sub_dir) -> str:
    return os.path.join(_rp_nmr_dir(), 'raw_pdbs', sub_dir)

def _rp_parsed_cifs_dir(sub_dir) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', sub_dir)


def rp_mmseqs_dir(sub_dir) -> str:
    return os.path.join(_rp_nmr_dir(), 'mmseqs', sub_dir)


def rp_mmseqs_fasta_dir(sub_dir) -> str:
    return os.path.join(rp_mmseqs_dir(sub_dir), 'fasta')


def rp_mmseqs_comb_fastas_dir(sub_dir) -> str:
    return os.path.join(rp_mmseqs_dir(sub_dir), 'combined_fastas')


def rp_mmseqs_results_dir(sub_dir) -> str:
    return os.path.join(rp_mmseqs_dir(sub_dir), 'results')


def rm_selfaligns_and_write_results(sub_dir: str, pdbid: str, pdf):
    """
    1. Remove self-alignments.
    2. Write non-empty mmseqs results pdf to csv.
    3. Write PDB id of empty results pdf to lst file.
    4. Return pdf which has had self-alignments removed.
    """
    pdf = pdf[pdf['query'] != pdf['target']]
    rp_mmseqs_dir_ = rp_mmseqs_dir(sub_dir)
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
            pdbids_no_idty = list(set(pdbids_no_idty))
            pdbids_no_idty.sort()
            with open(rp_zero_idty_lst, 'w') as f:
                f.write('\n'.join(pdbids_no_idty) + '\n')
        else:
            with open(rp_zero_idty_lst, 'w') as f:
                f.write(f'{pdbid}\n')
    return pdf


# RUNNING MMSEQS2 COMMANDS VIA SUBPROCESS.RUN():
def run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f: str):
    """
    Runs MMseqs2 all-vs-all search on the provided FASTA file.
    NOTE: NO NEED TO SPLIT THE PDBS INTO CHAINS AS INPUT, MMSEQS TAKES THE FASTA FILE OF EACH PDB WHICH CONTAINS
    ALL ITS CHAINS SEPARATED OUT AND RUNS EACH PDB-CHAIN AGAINST EVERY OTHER (INCLUDING ITSELF).

    Outputs:
      - results.m8 (tabular results)
      - Pandas DataFrame of alignments

    Note: Alignment results that fall below some internal threshold in MMseqs2 are just not returned at all.
    """
    fname = os.path.basename(rp_fasta_f).removesuffix('.fasta')
    rp_dir2delete = os.path.join(rp_mmseqs_dir(sub_dir=rp_fasta_f.split('/')[4]), 'dir2delete')
    shutil.rmtree(rp_dir2delete, ignore_errors=True)

    os.makedirs(rp_dir2delete, exist_ok=True)
    tmp_dir = os.path.join(rp_dir2delete, f'{fname}_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    output_m8 = os.path.join(rp_dir2delete, f'{fname}_results.m8')

    query_fasta = rp_fasta_f
    target_fasta = rp_fasta_f

    subprocess.run([
        'mmseqs', 'easy-search',
        query_fasta, target_fasta,     # all-vs-all
        output_m8, tmp_dir,
        '-s', '7.5',                   # high sensitivity
        '-e', '1000',                  # permissive e-value
        '--alignment-mode', '3',       # local alignment
        '--max-seqs', '10000',         # don't truncate results
        '--format-output', 'query,target,evalue,pident,alnlen,qcov,tcov'  # NOTE: ANY WHITE SPACES WILL BREAK THIS !
    ], check=True)

    pdf = pd.read_csv(output_m8, sep='\t', header=None,
                      names=['query', 'target', 'evalue', 'pident', 'alnlen', 'qcov', 'tcov'])
    return pdf  # shape=(60589, 7)


# def run_mmseqs2_createdb_cmd_all_vs_all(rp_fasta_f):
#     """
#     Runs MMseqs2 all-vs-all search on the provided FASTA file.
#     Outputs:
#       - results.m8 (tabular results)
#       - Pandas DataFrame of alignments
#     """
#     het_hom = rp_fasta_f.split('/')[4]
#     pdbid = os.path.basename(rp_fasta_f).removesuffix('.fasta')
#     rp_dir2delete = os.path.join(_rp_mmseqs_dir(het_hom), 'dir2delete')
#     shutil.rmtree(rp_dir2delete, ignore_errors=True)
#
#     os.makedirs(rp_dir2delete, exist_ok=True)
#     tmp_dir = os.path.join(rp_dir2delete, f'{pdbid}_tmp')
#     os.makedirs(tmp_dir, exist_ok=True)
#
#     db_name = os.path.join(rp_dir2delete, f'{pdbid}_db')
#     result_db = os.path.join(rp_dir2delete, f'{pdbid}_result')
#     output_m8 = os.path.join(rp_dir2delete, f'{pdbid}_results.m8')
#
#     # 1. Create database:
#     subprocess.run([
#         'mmseqs', 'createdb',
#         rp_fasta_f, db_name],
#         stdout=subprocess.DEVNULL, check=True)
#
#     # 2. Search: all-vs-all (with lower threshold):
#     subprocess.run([
#         'mmseqs', 'search',
#         db_name, db_name,
#         result_db, tmp_dir,
#         '--alignment-mode', '3',
#         '--threads', '4',
#         '-s', '7.5',
#         '-e', '1000'],
#         stdout=subprocess.DEVNULL, check=True)
#
#     # 3. Convert alignments to tabular output
#     subprocess.run([
#         'mmseqs', 'convertalis',
#         db_name, db_name,
#         result_db, output_m8,
#         '--format-output',
#         'query,target,evalue,pident,alnlen,qcov,tcov'],  # BEWARE SPACES BREAK THIS.
#         stdout=subprocess.DEVNULL, check=True)
#
#     return pd.read_csv(output_m8, sep=' ', header=None, names=['query, target, evalue, pident, alnlen, qcov, tcov'])
#

# def run_mmseqs2_across_all_het_or_hom_pdbids(het_hom):
#     rp_fasta_f = os.path.join(_rp_mmseqs_comb_fastas_dir(het_hom), 'het_combined_932_PDBids.fasta')
#     if het_hom == 'homomeric':
#         rp_fasta_f = os.path.join(_rp_mmseqs_comb_fastas_dir(het_hom), 'hom_combined_609_PDBids.fasta')
#     return run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f)


# VISUALISE FOR EMPIRICALLY DETERMINE EDITS TO THRESHOLDS:
def _view_pw_alignment_results(title: str, pdf_msq2):
    # Scatter: identity vs evalue (log scale)
    plt.figure(figsize=(10, 6))
    plt.scatter(pdf_msq2['pident'], pdf_msq2['evalue'], alpha=0.3)
    plt.yscale('log')
    plt.axvline(30, color='red', linestyle='--', label='30% identity threshold')
    plt.axhline(1e-3, color='green', linestyle='--', label='e-value = 1e-3')
    plt.xlabel('% Identity')
    plt.ylabel('E-value (log scale)')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def dedupe_rm_self_aligns(pdf):
    pdf = pdf[pdf['query'] != pdf['target']]  # (60589,7) --> (58008,7)
    # Create new col of pdbchain in query and target pair, to sort, to de-dupe:
    pdf['pair'] = pdf.apply(lambda row: tuple(sorted([row['query'], row['target']])), axis=1) # (58008,7) --> (58008,8)
    pdf = pdf.drop_duplicates(subset='pair').drop(columns='pair') # (58008,8) --> (33074,7)
    return pdf

def find_homologues_30_20_90(mmseq2_pdf, pid_chain: str):
    homologues_pdf = mmseq2_pdf[
        (mmseq2_pdf['evalue'] < 1e-3) &  # e-value threshold: < 0.001 (i.e. more significant matches)
        (mmseq2_pdf['pident'] >= 30.0) &  # sequence identity: >= 30%
        (mmseq2_pdf['alnlen'] >= 20) &  # alignment length: >= 20 residues
        (mmseq2_pdf['qcov'] >= 0.9) &  # query coverage: >=90%
        (mmseq2_pdf['tcov'] >= 0.9)  # target coverage: >= 90%
        ]
    rows_pidchain = homologues_pdf.loc[homologues_pdf['query'] == pid_chain]
    homologues = rows_pidchain['target'].tolist()
    rows_pidchain = homologues_pdf.loc[homologues_pdf['target'] == pid_chain]
    homologues = homologues + rows_pidchain['query'].tolist()
    return homologues


def calc_and_write_results_mmseqs2_all_vs_all(rp_fasta_f: str):
    """
    Note MMSeqs2 does not return results for alignments that are below a certain threshold.
    """
    pdf = run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f)
    print(f'pdf.shape={pdf.shape}')  # shape=(91333, 7)
    mmseqs2_result_pdf = pdf[pdf['query'] != pdf['target']]  # REMOVE ALL SELF-ALIGNMENTS
    print(f'mmseqs2_result_pdf.shape={mmseqs2_result_pdf.shape}')  # shape=(87944, 7)
    # Filter MMseqs2 hits based on the following homology thresholds that are heuristically chosen:
    homologues_pdf = pdf[
        (pdf['evalue'] < 1e-3) &    # e-value threshold: < 0.001 (i.e. more significant matches)
        (pdf['pident'] >= 30.0) &   # sequence identity: >= 30%
        (pdf['alnlen'] >= 20) &     # alignment length: >= 20 residues
        (pdf['qcov'] >= 0.9) &      # query coverage: >=90%
        (pdf['tcov'] >= 0.9)        # target coverage: >= 90%
        ]
    print(f'homologues_pdf.shape={homologues_pdf.shape}')
    _view_pw_alignment_results(title='homologues_30_20_90', pdf_msq2=homologues_pdf)
    non_homologues_pdf = pdf.drop(homologues_pdf.index)
    print(f'non_homologues_pdf.shape={non_homologues_pdf.shape}')
    _view_pw_alignment_results(title='non_homologues_30_20_90', pdf_msq2=non_homologues_pdf)

    dst_dir = rp_mmseqs_results_dir(sub_dir='hethom_combined')
    homologues_pdf.to_csv(os.path.join(dst_dir, 'homologues_30_20_90.csv'), index=False)  # shape=(39265, 7)
    non_homologues_pdf.to_csv(os.path.join(dst_dir, 'non_homologues_30_20_90.csv'), index=False)  # shape=(52068, 7)
    # The number of non_homologues is lower than the actual number because the alignments results returned
    # by MMSeqs2 excludes some that are just below some internal thresholds.
    return homologues_pdf, non_homologues_pdf


# FUNCTIONS FOR WRITING FASTA FILES:
def combine_all_het_and_hom_fastas_to_one_file():
    # destination dir and fasta file:
    rp_dst_hethom_comb_fasta_dir = rp_mmseqs_fasta_dir(sub_dir='hethom_combined')
    os.makedirs(rp_dst_hethom_comb_fasta_dir, exist_ok=True)
    rp_dst_combo_fastas_f = os.path.join(rp_dst_hethom_comb_fasta_dir, 'hethom_comb_1541_PDBids.fasta')

    # source fasta files:
    if not os.path.exists('../data/NMR/mmseqs/heteromeric/combined_fastas/het_all_932_PDBids.fasta'):
        print('het_all_932_PDBids.fasta not found')
    rp_het_fasta_f = os.path.join(rp_mmseqs_comb_fastas_dir('heteromeric'), 'het_all_932_PDBids.fasta')

    rp_hom_fasta_f = os.path.join(rp_mmseqs_comb_fastas_dir('homomeric'), 'hom_all_609_PDBids.fasta')

    with open(rp_dst_combo_fastas_f, 'wb') as dst:
        for src_path in [rp_het_fasta_f, rp_hom_fasta_f]:
            with open(src_path, 'rb') as src:
                shutil.copyfileobj(src, dst)


def combine_all_het_or_hom_fastas_into_one_file(het_hom: str):
    # destination dir and fasta file:
    rp_dst_comb_fastas_dir = rp_mmseqs_comb_fastas_dir(het_hom)
    os.makedirs(rp_dst_comb_fastas_dir, exist_ok=True)
    rp_dst_combo_fastas_f = os.path.join(rp_dst_comb_fastas_dir, 'het_combined_932_PDBids.fasta')

    if het_hom == 'homomeric':
        rp_dst_combo_fastas_f = os.path.join(rp_dst_comb_fastas_dir, 'hom_combined_609_PDBids.fasta')

    # source fasta files:
    rp_src_fasta_files = sorted(glob.glob(f'{rp_mmseqs_fasta_dir(het_hom)}/*.fasta'))

    with open(rp_dst_combo_fastas_f, 'w') as rp_dst_fasta_f:
        for rp_src_fasta_f in rp_src_fasta_files:
            with open(rp_src_fasta_f, 'r') as rp_src_f:
                for line in rp_src_f:
                    rp_dst_fasta_f.write(line)
    return rp_dst_combo_fastas_f


def write_fasta(het_hom: str, pdbid: str, chains: list) -> str:
    fasta_dir = rp_mmseqs_fasta_dir(het_hom)
    os.makedirs(fasta_dir, exist_ok=True)
    toks_dir = _rp_parsed_cifs_dir(het_hom)
    fasta_seq = []
    for chain in chains:
        pdbid_chain_pdf = pd.read_csv(os.path.join(toks_dir, f'{pdbid}_{chain}.ssv'), sep=' ')
        models = pdbid_chain_pdf['A_pdbx_PDB_model_num'].unique().tolist()
        pdf_one_model = pdbid_chain_pdf[pdbid_chain_pdf['A_pdbx_PDB_model_num'] == models[0]]
        aa3seq = pdf_one_model['S_mon_id'].tolist()
        aa1seq = map_aa3to1(aa3seq, pdbid, chain)
        fasta_seq.append(f'>{pdbid}_{chain}')
        fasta_seq.append(aa1seq)
    rp_fasta_file = os.path.join(fasta_dir, f'{pdbid}.fasta')

    with open(rp_fasta_file, 'w') as f:
        f.writelines('\n'.join(fasta_seq) + '\n')

    return rp_fasta_file



if __name__ == '__main__':
    rp_1A03 = '../data/NMR/mmseqs/homomeric/results/1A03.csv'
    pdf = pd.read_csv(rp_1A03)
    pdf = dedupe_rm_self_aligns(pdf)
    rp_fasta_f_ = os.path.join(rp_mmseqs_fasta_dir(sub_dir='hethom_combined'), 'hethom_comb_1541_PDBids.fasta')
    # result_pdf = run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f_)
    calc_and_write_results_mmseqs2_all_vs_all(rp_fasta_f_)

    # rp_pdbid_chains = os.path.join('..', 'data', 'NMR', 'multimodel_lists', 'het_multimod_2104_PidChains.txt')
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

    # def aa3to1_unused(pdbid: str, chain: str, three_char_seq: List[Tuple[int, str]]):
    #     aa_3to1 = {
    #         'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
    #         'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
    #         'TRP': 'W', 'TYR': 'Y'}
    #
    #     aa_3to1_natural_ptm = {
    #         'ASQ': 'S',  # Acetylserine
    #         'CGU': 'E',  # Gamma-carboxyglutamate
    #         'CIR': 'R',  # Citrulline
    #         'CSO': 'C',  # S-hydroxycysteine
    #         'CSX': 'C',  # Cysteine S-oxide
    #         'FME': 'M',  # N-formylmethionine
    #         'HIC': 'H',  # 4-methylhistidine
    #         'HYP': 'P',  # Hydroxyproline
    #         'KCX': 'K',  # N-epsilon-carboxylysine
    #         'MLY': 'K',  # Alternate for methyllysine
    #         'MLZ': 'K',  # N6-methyllysine
    #         'MSE': 'M',  # Selenomethionine
    #         'PCA': 'Q',  # Pyroglutamic acid
    #         'PTR': 'Y',  # Phosphotyrosine
    #         'PYR': 'Q',  # Alternate code for pyroglutamic acid
    #         'SEP': 'S',  # Phosphoserine
    #         'TPO': 'T',  # Phosphothreonine
    #     }
    #
    #     aa_3to1_manmade = {
    #         'DAL': 'A',  # D-alanine
    #         'DAR': 'R',  # D-arginine
    #         'DCY': 'C',  # D-cysteine
    #         'DGL': 'E',  # D-glutamic acid
    #         'DGN': 'Q',  # D-glutamine
    #         'DHI': 'H',  # D-histidine
    #         'DIL': 'I',  # D-isoleucine
    #         'DLE': 'L',  # D-leucine
    #         'DLY': 'K',  # D-lysine
    #         'DPH': 'F',  # D-phenylalanine
    #         'DSN': 'S',  # D-serine
    #         'DTH': 'T',  # D-threonine
    #         'DTY': 'Y',  # D-tyrosine
    #         'DVA': 'V',  # D-valine
    #     }
    #
    #     fasta_seq = []
    #     for i, (_, res) in enumerate(three_char_seq):
    #         if res not in aa_3to1:
    #             # raise KeyError(f'{res} not found in list of three-letter codes. PDBid={pdbid}, Chain={cif_chain}')
    #             print(f'{pdbid}_{chain}: {res} not one of the 20 natural residues.')
    #             if res not in aa_3to1_natural_ptm:
    #                 print(f'{pdbid}_{chain}: {res} also not in my list of 17 post-translationally modified residues.')
    #                 if res not in aa_3to1_manmade:  # not natural and not PTM - so remove:
    #                     print(f'{pdbid}_{chain}: {res} also not in my list of 14 synthetic residues. '
    #                           f'Therefore this pdbid_chain will be removed, as I have no idea what it is !')
    #                     return None
    #                 else:
    #                     print(f'{pdbid}_{chain}: {res} is in my list of 14 synthetic residues.')
    #                     print(f'Therefore this pdbid_chain will be removed, as it contains synthetic residues.')
    #                     return None
    #             else:  # natural PTM residue:
    #                 aa1char = aa_3to1_natural_ptm[res]
    #                 print(f'{pdbid}_{chain}: {res} is in PTM residues list and is being mapped to {aa1char}')
    #                 if (i + 1) % 80 == 0:
    #                     fasta_seq.append(aa1char + '\n')
    #                 else:
    #                     fasta_seq.append(aa1char)
    #         else:  # natural residue:
    #             if (i + 1) % 80 == 0:
    #                 fasta_seq.append(aa_3to1[res] + '\n')
    #             else:
    #                 fasta_seq.append(aa_3to1[res])
    #
    #     return ''.join(fasta_seq)
