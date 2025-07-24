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

def _rp_nmr_dir():
    return os.path.join('..', 'data', 'NMR')


def _rp_rawcifs_dir(het_hom: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'raw_cifs', het_hom)


def _rp_tokcifs_dir(het_hom) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', het_hom)


def _rp_mmseqs_dir(het_hom) -> str:
    return os.path.join(_rp_nmr_dir(), 'mmseqs', het_hom)


def _rp_mmseqs_fasta_dir(het_hom) -> str:
    return os.path.join(_rp_mmseqs_dir(het_hom), 'fasta')


def _rp_mmseqs_comb_fastas_dir(het_hom) -> str:
    return os.path.join(_rp_mmseqs_dir(het_hom), 'combined_fastas')


def _rp_mmseqs_hethom_comb_fasta_dir() -> str:
    return os.path.join(_rp_nmr_dir(), 'mmseqs', 'hethom', 'fasta')


def combine_het_hom_fastas_to_1_file():
    # destination dir and fasta file:
    rp_dst_hethom_comb_fasta_dir = _rp_mmseqs_hethom_comb_fasta_dir()
    os.makedirs(rp_dst_hethom_comb_fasta_dir, exist_ok=True)
    rp_dst_combo_fastas_f = os.path.join(rp_dst_hethom_comb_fasta_dir, 'hethom_comb_1541_PDBids.fasta')

    # source fasta files:
    if not os.path.exists('../data/NMR/mmseqs/heteromeric/combined_fastas/het_all_932_PDBids.fasta'):
        print('het_all_932_PDBids.fasta not found')
    rp_het_fasta_f = os.path.join(_rp_mmseqs_comb_fastas_dir('heteromeric'), 'het_all_932_PDBids.fasta')

    rp_hom_fasta_f = os.path.join(_rp_mmseqs_comb_fastas_dir('homomeric'), 'hom_all_609_PDBids.fasta')

    with open(rp_dst_combo_fastas_f, 'wb') as dst:
        for src_path in [rp_het_fasta_f, rp_hom_fasta_f]:
            with open(src_path, 'rb') as src:
                shutil.copyfileobj(src, dst)


def combine_all_fastas_in_1_file(het_hom: str):
    # destination dir and fasta file:
    rp_dst_comb_fastas_dir = _rp_mmseqs_comb_fastas_dir(het_hom)
    os.makedirs(rp_dst_comb_fastas_dir, exist_ok=True)
    rp_dst_combo_fastas_f = os.path.join(rp_dst_comb_fastas_dir, 'het_combined_932_PDBids.fasta')
    if het_hom == 'homomeric':
        rp_dst_combo_fastas_f = os.path.join(rp_dst_comb_fastas_dir, 'hom_combined_609_PDBids.fasta')

    # source fasta files:
    rp_src_fasta_files = sorted(glob.glob(f'{_rp_mmseqs_fasta_dir(het_hom)}/*.fasta'))

    with open(rp_dst_combo_fastas_f, 'w') as rp_dst_fasta_f:
        for rp_src_fasta_f in rp_src_fasta_files:
            with open(rp_src_fasta_f, 'r') as rp_src_f:
                for line in rp_src_f:
                    rp_dst_fasta_f.write(line)

    return rp_dst_combo_fastas_f


def write_fasta(het_hom: str, pdbid: str, chains: list) -> str:
    fasta_dir = _rp_mmseqs_fasta_dir(het_hom)
    os.makedirs(fasta_dir, exist_ok=True)
    toks_dir = _rp_tokcifs_dir(het_hom)
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


def run_easysearch_mmseqs2_allvsall(rp_fasta_f):
    """
    Runs MMseqs2 all-vs-all search on the provided FASTA file.
    Outputs:
      - results.m8 (tabular results)
      - Pandas DataFrame of alignments
    """
    het_hom = rp_fasta_f.split('/')[4]
    pdbid = os.path.basename(rp_fasta_f).removesuffix('.fasta')
    rp_dir2delete = os.path.join(_rp_mmseqs_dir(het_hom), 'dir2delete')
    shutil.rmtree(rp_dir2delete, ignore_errors=True)

    os.makedirs(rp_dir2delete, exist_ok=True)
    tmp_dir = os.path.join(rp_dir2delete, f'{pdbid}_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    output_m8 = os.path.join(rp_dir2delete, f'{pdbid}_results.m8')

    query_fasta = rp_fasta_f
    target_fasta = rp_fasta_f

    subprocess.run([
        'mmseqs', 'easy-search',
        query_fasta, target_fasta,         # all-vs-all
        output_m8, tmp_dir,
        '-s', '7.5',                    # high sensitivity
        '-e', '1000',                   # permissive e-value
        '--alignment-mode', '3',       # local alignment
        '--max-seqs', '10000',         # don't truncate results
        '--format-output', 'query,target,evalue,pident,alnlen'
    ], check=True)

    return pd.read_csv(output_m8, sep='\t', header=None, names=['query', 'target', 'evalue', 'pident', 'alnlen'])

def run_mmseqs2_all_vs_all(rp_fasta_f):
    """
    Runs MMseqs2 all-vs-all search on the provided FASTA file.
    Outputs:
      - results.m8 (tabular results)
      - Pandas DataFrame of alignments
    """
    het_hom = rp_fasta_f.split('/')[4]
    pdbid = os.path.basename(rp_fasta_f).removesuffix('.fasta')
    rp_dir2delete = os.path.join(_rp_mmseqs_dir(het_hom), 'dir2delete')
    shutil.rmtree(rp_dir2delete, ignore_errors=True)

    os.makedirs(rp_dir2delete, exist_ok=True)
    tmp_dir = os.path.join(rp_dir2delete, f'{pdbid}_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    db_name = os.path.join(rp_dir2delete, f'{pdbid}_db')
    result_db = os.path.join(rp_dir2delete, f'{pdbid}_result')
    output_m8 = os.path.join(rp_dir2delete, f'{pdbid}_results.m8')

    # 1. Create database:
    subprocess.run([
        'mmseqs', 'createdb',
        rp_fasta_f, db_name],
        stdout=subprocess.DEVNULL, check=True)

    # 2. Search: all-vs-all (with lower threshold):
    subprocess.run([
        'mmseqs', 'search',
        db_name, db_name,
        result_db, tmp_dir,
        '--alignment-mode', '3',
        '--threads', '4',
        '-s', '7.5',
        '-e', '1000'],
        stdout=subprocess.DEVNULL, check=True)

    # 3. Convert alignments to tabular output
    subprocess.run([
        'mmseqs', 'convertalis',
        db_name, db_name,
        result_db, output_m8,
        '--format-output',
        'query,target,evalue,pident,alnlen'],
        stdout=subprocess.DEVNULL, check=True)

    return pd.read_csv(output_m8, sep=' ', header=None, names=['query','target','evalue','pident','alnlen'])


def filter_and_write_results(het_hom: str, pdbid: str, pdf):
    pdf = pdf[pdf['query'] != pdf['target']]
    rp_mmseqs_dir_ = _rp_mmseqs_dir(het_hom)
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


def filter_out_100_pident(pdf):
    return pdf[pdf['pident'] < 100]

def run_mmseqs2_across_all_pdbids(het_hom):
    rp_fasta_f = os.path.join(_rp_mmseqs_comb_fastas_dir(het_hom), 'het_combined_932_PDBids.fasta')
    if het_hom == 'homomeric':
        rp_fasta_f = os.path.join(_rp_mmseqs_comb_fastas_dir(het_hom), 'hom_combined_609_PDBids.fasta')
    return run_easysearch_mmseqs2_allvsall(rp_fasta_f)

def run_mmseqs2_across_all_pdbids_combined():
    rp_fasta_f = os.path.join(_rp_mmseqs_hethom_comb_fasta_dir(), 'hethom_comb_1541_PDBids.fasta')
    return run_easysearch_mmseqs2_allvsall(rp_fasta_f)


if __name__ == '__main__':

    run_mmseqs2_across_all_pdbids_combined()
    pass

    # combine_het_hom_fastas_to_1_file()
    # het_hom = 'heteromeric'
    # txt_f = 'het_multimod_2104_PidChains.txt'
    # fasta_f = 'het_all_932_PDBids.fasta'

    # het_hom = 'homomeric'
    # txt_f = 'hom_multimod_1421_PidChains.txt'
    # fasta_f = 'hom_all_609_PDBids.fasta'
    # pdbids, pdbid_chains_dict = pms._read_multimodel_pdbid_chains(txt_f)

    # for pid, chains in pdbid_chains_dict.items():
    #     # write_fasta(het_hom=het_hom, pdbid=pid, chains=chains)
    #     rp_fasta_file = os.path.join(_rp_mmseqs_fasta_dir(het_hom), f'{pid}.fasta')
    #     result_pdf = run_easysearch_mmseqs2_allvsall(rp_fasta_file)
    #     pdf_ = filter_and_write_results(het_hom, pdbid=pid, pdf=result_pdf)




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
