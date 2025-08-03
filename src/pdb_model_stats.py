import os, glob
from time import time
from typing import Tuple
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio import SeqIO
import RMSD
import mmseqs2


def _total_chain_count_and_year(het_hom: str, pdbid: str, use_mmcif: bool) -> tuple:
    """
    # cif_dict = MMCIF2Dict.MMCIF2Dict(rpath_cif)
    # year = cif_dict.get('_citation.year', ['NA'])[0]
    # A few of these return "?", so I chose to use 'MMCIFParser().get_structure(cif)' instead,
    # which often gives same year but sometimes differs by as much as ~3 years.
    """
    parser = MMCIFParser(QUIET=True)  # 'QUIET' is for debugging possible oddities that might pop up in PDB records.
    if not use_mmcif:
        parser = PDBParser(QUIET=True)
    rp_raw_cif_dir = os.path.join('..', 'data', 'NMR', 'raw_cifs', het_hom)
    rp_cif = os.path.join(rp_raw_cif_dir, f'{pdbid}.cif')
    bio_struct = parser.get_structure('', rp_cif)
    total_chain_count = len(list(next(bio_struct.get_models()).get_chains()))
    year = bio_struct.header['deposition_date'][:4]
    return total_chain_count, year


def _pdbid_dict_chain(rp_pdbid_chains: str) -> Tuple[list, dict]:
    """
    Takes relative path of txt file that has list of sol NMR PDBid_chains with > 1 model.
    Makes a dict with unique PDBid as the key, mapped to all its chains, according to the list in the txt file.
    Returns a list of the PDBid-chains, and a dictionary mapping PDB ids to chains.
    """
    # READ IN THE 2104 HETEROMERIC OR 1421 HOMOMERIC PDBid_chains.
    with open(rp_pdbid_chains, 'r') as f:
        # pdbid_chains = [line.rstrip('\n') for line in f]
        pdbid_chains = f.read().split()

    # MAKE DICT OF THE 932 HETEROMERIC OR 609 HOMOMERIC PDB ids, MAPPED TO THEIR CORRESPONDING CHAINS.
    pdbid_chains_dict = defaultdict(list)
    pdbid_chains.sort()  # expecting already sorted (i.e. before writing out), but to avoid relying on other code.
    for pdbid_chain in pdbid_chains:
        pid, ch = pdbid_chain.split('_')
        pdbid_chains_dict[pid].append(ch)
    return pdbid_chains, dict(pdbid_chains_dict)


def _read_multimodel_pdbid_chains(txt_f: str) -> Tuple[list, dict, list]:
    rp_pdbid_chains = os.path.join('..', 'data', 'NMR', 'multimodel_lists', txt_f)
    pidChains_list, pidChains_dict = _pdbid_dict_chain(rp_pdbid_chains)
    pdbids = list(pidChains_dict.keys())
    return pdbids, pidChains_dict, pidChains_list




def _calc_identity_for_stats(het_hom: str, pdbid: str, filt_pdf) -> str:
    if not filt_pdf.empty:
        print(f'There is some heterogeneity of sequence between some chains. '
              f'See data/NMR/mmseqs/{het_hom}/results/{pdbid}.csv')
        identity = '(chains not identical, see mmseqs pdf)'
    else:
        print('All chains have high identity, if not 100%')
        identity = '0'
    return identity


def _calc_mmseqs2_single_pdb(het_hom: str, pid: str):
    rp_fasta_f = os.path.join(mmseqs2.rp_mmseqs_fasta_dir(het_hom), f'{pid}.fasta')
    result_pdf = mmseqs2.run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f)
    return result_pdf


def generate_stats(het_hom: str, use_mmcif=True):
    """
    After parsing raw mmCIFs, the resulting PDBid_chains were written to a .lst file, which is used as follows:

    1. Read raw mmCIFs via Bio.PDB.MMCIFParser for each PDBid in .lst to get:
        - deposition year for each PDBid
        - total number of chains for each PDBid.

    2. Read PDBid_chains .lst file to get:
        - number of protein chains per PDBid.

    3. Read parsed SSV files to dataframes for each PDBid_chain in .lst to get:
        - number of models per PDBid_chain.
        - number of alpha carbons per PDBid_chain.

    4. Calculate sequence identity by MMseqs2 between:

        - every chain of each PDBid from list.
        - every PDBid_chain with every other.

    5. Calculate RMSD for each model against the average of those models.
    """
    start = time()
    stats = []

    if het_hom == 'heteromeric':
        pdbids, pidChains_dict, pidChains_list = _read_multimodel_pdbid_chains('het_multimod_2104_PidChains.txt')
    if het_hom == 'homomeric':
        pdbids, pidChains_dict, pidChains_list = _read_multimodel_pdbid_chains('hom_multimod_1421_PidChains.txt')

    rp_parsed_cif_dir = os.path.join('..', 'data', 'NMR', 'parsed_cifs', het_hom)

    for pid, chains in pidChains_dict.items():
        total_chain_count, year = _total_chain_count_and_year(het_hom, pid, use_mmcif)
        mmseqs_pdf = _calc_mmseqs2_single_pdb(het_hom, pid)

        for chain in chains:
            pid_chain = f'{pid}_{chain}'
            identity = _calc_identity_for_stats(het_hom, pid, mmseqs_pdf)
            rmsds, model_nums, one_pdbidchain_in_list = RMSD.rmsds_across_models_in_one_pidChain(het_hom, pid_chain, rp_parsed_cif_dir)
            max_rmsd, min_rmsd = np.max(rmsds), np.min(rmsds)  # this line is just for debugging
            rmsd_pdf = pd.DataFrame({'pdbid_chain': one_pdbidchain_in_list, 'model_num':model_nums, 'rmsd': rmsds})
            rp_rmsd_dir = os.path.join('..', 'data', 'NMR', 'RMSD', het_hom)
            os.makedirs(rp_rmsd_dir, exist_ok=True)
            rp_rmsd_csv = os.path.join(rp_rmsd_dir, f'{pid_chain}.csv')
            rmsd_pdf.to_csv(rp_rmsd_csv, index=False)
            rp_tok_cif = os.path.join(rp_parsed_cif_dir, f'{pid_chain}.ssv')
            pdbid_chain_pdf = pd.read_csv(rp_tok_cif, sep=' ')
            model_count = len(pdbid_chain_pdf['A_pdbx_PDB_model_num'].unique())
            ca_count = pdbid_chain_pdf.shape[0] / model_count

            stats.append({
                'PDBid': pid,
                'het_hom': het_hom,
                'year': year,
                '#model': model_count,
                '#chains': total_chain_count,
                '#protChains': len(pidChains_dict[pid]),
                'chain': chain,
                '#alphaCarbs': ca_count,
                'identity': identity,
                'minRMSD': np.min(rmsds),
                'maxRMSD': np.max(rmsds),
                # 'similarity': 100  # TODO
            })
    pdf = pd.DataFrame(stats)

    pdf_sorted = pdf.sort_values(by=['year', '#model'], ascending=[True, True])
    print(pdf_sorted.dtypes)
    for col_to_cast in ['PDBid', 'chain']:
        pdf_sorted[col_to_cast] = pdf_sorted[col_to_cast].astype('string')
    print(pdf_sorted.dtypes)
    protein_lengths(pdf_sorted[['#alphaCarbs']].values)
    num_cifs = int(pidChains_list.split('_')[2])
    print(f'Completed {num_cifs} CIFs in {round((time() - start) / 60)} mins')
    return pdf_sorted


def protein_lengths(ca_count):
    # Example data: replace this with your actual ca_count data
    # ca_count = np.random.randint(50, 500, size=100)
    sorted_lengths = np.sort(ca_count)
    x = np.arange(len(sorted_lengths))
    plt.figure(figsize=(10, 6))
    plt.plot(x, sorted_lengths, marker='o', linestyle='-', linewidth=1)
    plt.xlabel('PDBs')
    plt.ylabel('Protein length (C-alpha count)')
    plt.title('Distribution of Protein Lengths (Sorted)')
    plt.xticks([])
    plt.tight_layout()
    plt.show()


def fasta_size_distribution(rp_fasta_f, x_limit_200=False):
    seq_lengths = [len(record.seq) for record in SeqIO.parse(rp_fasta_f, 'fasta')]
    length_counts = Counter(seq_lengths)

    sorted_lengths = sorted(length_counts.keys())
    counts = [length_counts[length] for length in sorted_lengths]

    fig, ax1 = plt.subplots(figsize=(10, 5))

    # Plot exact counts
    ax1.plot(sorted_lengths, counts, color='black', linewidth=1.5, label='Exact Count')
    ax1.set_xlabel('Sequence length', fontsize=10)
    ax1.set_ylabel('Number of sequences', fontsize=10)
    ax1.tick_params(axis='y', labelsize=9)

    # Aesthetic adjustments for primary axis
    ax1.set_xticks(np.arange(0, max(sorted_lengths) + 1, 20))
    ax1.set_xticklabels(np.arange(0, max(sorted_lengths) + 1, 20), rotation=-90, fontsize=8)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(True, linestyle=':', linewidth=0.5, color='lightgray')

    if x_limit_200:
        ax1.set_xlim(0, 200)

    # Secondary axis for KDE
    ax2 = ax1.twinx()
    sns.kdeplot(
        seq_lengths,
        ax=ax2,
        color='blue',
        linestyle='--',
        linewidth=1.5,
        label='KDE',
        bw_adjust=0.4,
        cut=0
    )
    ax2.set_ylabel('Density (KDE)', fontsize=10)
    ax2.tick_params(axis='y', labelsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Title and layout
    plt.title('Distribution of PDB Sequence Lengths', fontsize=12)
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':

    # _calc_mmseqs2(het_hom='heteromeric', pid='0_all_het', chains=[])
    # _calc_mmseqs2(het_hom='homomeric', pid='0_all_hom', chains=[])
    # _calc_mmseqs2(het_hom='', pid='', chains=[])

    # _calc_mmseqs2(het_hom='heteromeric', pdbid='1AI0',
    #               chains=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'])
    rp_fasta_f_ = os.path.join(mmseqs2.rp_mmseqs_fasta_dir(het_hom='hethom_combined'),
                               'hethom_comb_1541_PDBids.fasta')
    fasta_size_distribution(rp_fasta_f_, x_limit_200=True)

    # rp_het_fasta_f_ = os.path.join(mmseqs2.rp_mmseqs_comb_fastas_dir(het_hom='heteromeric'), 'het_all_932_PDBids.fasta')
    # rp_hom_fasta_f_ = os.path.join(mmseqs2.rp_mmseqs_comb_fastas_dir(het_hom='homomeric'), 'hom_all_609_PDBids.fasta')
    # stats_pdf = generate_stats(rp_het_fasta_f_, rp_hom_fasta_f_)
    # pass


    # # Columns to blank out when duplicated:
    # cols_to_blank = ['PDBid', 'het_hom', 'year', '#model', '#chain', '#prot_chain', 'identity', 'similarity']
    # # For each column, within each PDBid group, replace duplicates with ""
    # for col in cols_to_blank:
    #     pdf[col] = pdf.groupby('PDBid')[col].transform(lambda x: x.mask(x.duplicated()))
    #
    # cols_to_cast_int = ['year', '#chain', '#model', '#prot_chain']
    # for col in cols_to_cast_int:
    #     pdf[col] = pdf[col].astype('Int64')
    # print(pdf)
    # # pdf = pdf.fillna('')
    # csv_dst_dir = f'../data/NMR/stats/{_meric}'
    # os.makedirs(csv_dst_dir, exist_ok=True)
    # pdf.to_csv(f'{csv_dst_dir}/from_{len(cifs)}_cifs.csv', index=False)

    # TODO generate 2d plots (akin to bar graphs) of datasets for:
    #   - protein lengths (CA count)
    #   - rmsd min and max values
    #   - number of models
    #   - date of depoosition

