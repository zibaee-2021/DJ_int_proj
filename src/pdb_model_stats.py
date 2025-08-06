import os, glob
from time import time
from typing import Tuple
import math
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

# BUILDING RELATIVE PATHS:
def _rp_nmr_dir():
    return os.path.join('..', 'data', 'NMR')

def _rp_mmseqs_dir(sub_dir) -> str:
    return os.path.join(_rp_nmr_dir(), 'mmseqs', sub_dir)

def _rp_mmseqs_fasta_dir(sub_dir) -> str:
    return os.path.join(_rp_mmseqs_dir(sub_dir), 'fasta')

def _total_chain_count_and_year(rp_raw_struct_f: str, biopython_parser) -> Tuple[int, int]:
    """
    # The following code here in the docstring is not used, for the reason given below:
    # cif_dict = MMCIF2Dict.MMCIF2Dict(rpath_cif)
    # year = cif_dict.get('_citation.year', ['NA'])[0]
    # A few of these return "?", so I chose to use 'MMCIFParser().get_structure(cif)' instead,
    # which often gives same year but sometimes differs by as much as ~3 years.
    """
    # Note: Parser.get_structure() is quite slow, 15-20 secs:
    bio_struct = biopython_parser.get_structure('', rp_raw_struct_f)
    total_chain_count = len(list(next(bio_struct.get_models()).get_chains()))
    year = bio_struct.header['deposition_date'][:4]
    return total_chain_count, int(year)


def _pdbid_dict_chain(pdbid_chains: list) -> Tuple[list, dict]:
    """
    Takes relative path of txt file that has list of sol NMR PDBid_chains with > 1 model.
    Makes a dict with unique PDBid as the key, mapped to all its chains, according to the list in the txt file.
    Returns a list of the PDBid-chains, and a dictionary mapping PDB ids to chains.
    """
    pdbid_chains_dict = defaultdict(list)
    pdbid_chains.sort()
    for pdbid_chain in pdbid_chains:
        pid, ch = pdbid_chain.split('_')
        pdbid_chains_dict[pid].append(ch)
    return pdbid_chains, dict(pdbid_chains_dict)


def _read_multimodel_pdbid_chains(rp_pidChains: str) -> Tuple[list, dict, list]:

    if rp_pidChains.endswith('.txt'):
        with open(rp_pidChains, 'r') as f:
            # pdbid_chains = [line.rstrip('\n') for line in f]
            pdbid_chains = f.read().split()
    else:  # .lst file
        with open(rp_pidChains, 'r') as f:
            pdbid_chains = f.read().splitlines()

    pdbid_chains.sort()
    pidChains_list, pidChains_dict = _pdbid_dict_chain(pdbid_chains)
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


def generate_stats(sub_dir: str, rp_pidchains_lst_f: str, rp_fasta_f: str, run_and_write_mmseqs2=False,
                   run_and_write_rmsd=False, use_mmcif=True):
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

    6. Calculate TM-score for homologous PDBchains. (Although, all vs all might be done and included anyway).
    """
    start = time()
    stats = []
    mmseqs_pdf = None
    fname = os.path.basename(rp_fasta_f).removesuffix('.fasta')
    rp_mmseqs2_results_dir = os.path.join(_rp_mmseqs_dir(sub_dir), 'results')
    os.makedirs(rp_mmseqs2_results_dir, exist_ok=True)
    rp_mmseqs2_csv_f = os.path.join(rp_mmseqs2_results_dir, f'{fname}.csv')

    if use_mmcif:
        pdbcif = 'cif'
        parser = MMCIFParser(QUIET=True)  # Used in _total_chain_count_and_year() below (~line 150).
        # 'QUIET' is for debugging any problems that might pop up in PDB records.
    else:
        pdbcif = 'pdb'
        parser = PDBParser(QUIET=True)
    rp_raw_struct_het_dir = os.path.join(_rp_nmr_dir(), f'raw_{pdbcif}s', 'heteromeric')
    rp_raw_struct_hom_dir = os.path.join(_rp_nmr_dir(), f'raw_{pdbcif}s', 'homomeric')
    rp_raw_het_struct_files = glob.glob(os.path.join(rp_raw_struct_het_dir, f'*.{pdbcif}'))
    rp_raw_hom_struct_files = glob.glob(os.path.join(rp_raw_struct_hom_dir, f'*.{pdbcif}'))
    rp_raw_struct_files = rp_raw_het_struct_files + rp_raw_hom_struct_files
    raw_pids = [os.path.basename(rp_raw_struct_f).removesuffix(f'.{pdbcif}') for rp_raw_struct_f in rp_raw_struct_files]
    raw_het_pids = [os.path.basename(rp_raw_het_struct_f).removesuffix(f'.{pdbcif}') for rp_raw_het_struct_f in rp_raw_het_struct_files]

    if run_and_write_mmseqs2:
        mmseq2_pdf = mmseqs2.run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f)
        mmseq2_pdf = mmseqs2.dedupe_rm_self_aligns(mmseq2_pdf)
        mmseq2_pdf.to_csv(rp_mmseqs2_csv_f, index=False)  # shape=(33074, 7)
    else:  # This relies on it having been run, filtered and written out separately beforehand, so can just read in:
        mmseq2_pdf = pd.read_csv(rp_mmseqs2_csv_f)  # shape=(33074, 7)

    pdbids, pidchains_dict, pidchains_list = _read_multimodel_pdbid_chains(rp_pidchains_lst_f)

    rp_parsed_cif_dir = os.path.join(_rp_nmr_dir(), 'parsed_cifs', sub_dir)
    rp_raw_struct_dir = os.path.join(_rp_nmr_dir(), f'raw_{pdbcif}s', 'hethom_combined')

    rp_rmsd_mean_coords_dir = os.path.join(_rp_nmr_dir(), 'RMSD', sub_dir, 'mean_coords')

    for i, (pid, chains) in enumerate(pidchains_dict.items()):
        if i < 1507:
            continue

        rp_raw_struct_f = os.path.join(rp_raw_struct_dir, f'{pid}.{pdbcif}')

        total_chain_count, year = _total_chain_count_and_year(rp_raw_struct_f, parser)

        for chain in chains:
            pid_chain = f'{pid}_{chain}'
            rp_parsed_cifs_ssv = os.path.join(rp_parsed_cif_dir, f'{pid}_{chain}.ssv')
            rp_mean_coords_csv = os.path.join(rp_rmsd_mean_coords_dir, f'{pid}_{chain}.csv')
            if run_and_write_rmsd:
                rmsds, model_nums, pdc4pdf = RMSD.calc_rmsds_of_models(rp_parsed_cifs_ssv, rp_mean_coords_csv)
                rmsd_pdf = pd.DataFrame({'pdbid_chain': pdc4pdf, 'model_num': model_nums, 'rmsd': rmsds})
                rp_rmsd_dir = os.path.join(_rp_nmr_dir(), 'RMSD', sub_dir)
                os.makedirs(rp_rmsd_dir, exist_ok=True)
                rp_rmsd_csv = os.path.join(rp_rmsd_dir, f'{pid_chain}.csv')
                rmsd_pdf.to_csv(rp_rmsd_csv, index=False)
            else: # This relies on it having been run separately beforehand, and we can just read it in now:
                rp_rmsd_per_model_dir = os.path.join(_rp_nmr_dir(), 'RMSD', sub_dir)
                rp_rmsd_per_model_csv_f = os.path.join(rp_rmsd_per_model_dir, f'{pid_chain}.csv')
                rmsd_pdf = pd.read_csv(rp_rmsd_per_model_csv_f)
            rmsds = rmsd_pdf[['rmsd']].values
            if len(rmsds) > 0:
                min_rmsd, max_rmsd = np.min(rmsds), np.max(rmsds)
            else:
                min_rmsd, max_rmsd = np.nan, np.nan
            rp_parsed_cif_ssv_f = os.path.join(rp_parsed_cif_dir, f'{pid_chain}.ssv')
            pdbid_chain_pdf = pd.read_csv(rp_parsed_cif_ssv_f, sep=' ')
            model_count = len(pdbid_chain_pdf['A_pdbx_PDB_model_num'].unique())
            ca_count = int(pdbid_chain_pdf.shape[0] / model_count)
            homologues = mmseqs2.find_homologues_30_20_90(mmseq2_pdf, pid_chain)

            if use_mmcif:
                assert pid in raw_pids, f'{pid}.cif not round in raw_cifs.'
            else:
                if pid not in ['9D9B', '7ZE0', '9D9C', '9D9A']: # (These 4 are not found in 'legacy' PDB files on RCSB)
                    assert pid in raw_pids, f'{pid}.pdb not round in raw_pdbs.'

            het_hom = 'het' if pid in raw_het_pids else 'hom'
            stats.append({
                'PDBid': pid,
                'het_hom': het_hom,
                'year': year,
                '#model': model_count,
                '#chains': total_chain_count,
                '#protChains': len(pidchains_dict[pid]),
                'chain': chain,
                '#alphaCarbs': ca_count,
                'homologues': homologues,
                'minRMSD': min_rmsd,
                'maxRMSD': max_rmsd,
                # 'similarity': 100  # TODO
            })
    pdf = pd.DataFrame(stats)

    pdf_sorted = pdf.sort_values(by=['year', '#model'], ascending=[True, True])
    print(pdf_sorted.dtypes)
    for col_to_cast in ['PDBid', 'chain']:
        pdf_sorted[col_to_cast] = pdf_sorted[col_to_cast].astype('string')
    print(pdf_sorted.dtypes)
    num_pdbs = len(pidchains_list)
    print(f'Completed {num_pdbs} {pdbcif}s in {round((time() - start) / 60)} mins')
    return pdf_sorted


def plot_protein_lengths_binned(ca_counts: list, bin_size: int = 1):
    ca_counts = np.array(ca_counts)

    if bin_size == 1:
        heights = ca_counts
        x = np.arange(len(ca_counts))
    else:
        n_bins = len(ca_counts) // bin_size
        trimmed = ca_counts[:n_bins * bin_size]
        grouped = trimmed.reshape(n_bins, bin_size)
        heights = grouped.mean(axis=1)
        x = np.arange(n_bins)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x, heights, color='lightgrey', edgecolor='gainsboro', linewidth=0.5, width=1.0, align='edge')
    ax.set_xlim(left=-5)

    ax.set_xlabel('PDBchains' if bin_size == 1 else f'Binned PDBchains (bin size = {bin_size})', fontsize=10)
    ax.set_xticklabels(np.arange(0, len(ca_counts) + 1, 20), rotation=-90, fontsize=8)
    ax.set_ylabel('CA count', fontsize=10)
    ax.set_title('Number of CA in each PDBchain', fontsize=12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('lightgrey')
    ax.spines['bottom'].set_color('lightgrey')

    target_n_ticks = 30
    tick_interval = max(1, int(math.ceil(len(x) / target_n_ticks / 10.0)) * 10)
    ax.set_xticks(np.arange(0, len(x) + 1, tick_interval))

    fig.tight_layout()
    plt.show()


def _calc_ca_counts(rp_parsed_cifs_ssvs: list) -> list:
    ca_counts = list()
    for rp_parsed_cif_ssv in rp_parsed_cifs_ssvs:
        pdf = pd.read_csv(rp_parsed_cif_ssv, sep=' ')
        ca_counts_all_models = pdf.shape[0]
        model_count = len(pdf['A_pdbx_PDB_model_num'].unique())
        ca_count = int(ca_counts_all_models / model_count)
        if ca_count <= 3:
            print(f'Less than 4 CAs! {os.path.basename(rp_parsed_cif_ssv).removesuffix('.ssv')}={ca_count}')
        ca_counts.append(ca_count)
    return sorted(ca_counts)


def plot_fasta_size_distribution(rp_fasta_f, x_limit_220=False):
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

    if x_limit_220:
        ax1.set_xlim(0, 220)
    else:
        ax1.set_xlim(0, max(sorted_lengths) + 1)

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
    # VISUALISATIONS:
    rp_fasta_f_ = os.path.join(_rp_mmseqs_fasta_dir(sub_dir='multimod_2713_hetallchains_hom1chain'),
                               'multimod_2713_hetallchains_hom1chain.fasta')
    # plot_fasta_size_distribution(rp_fasta_f_, x_limit_220=False)
    rp_parsed_cifs_ssvs_ = sorted(glob.glob(os.path.join(_rp_nmr_dir(), 'parsed_cifs',
                                                 'multimod_2713_hetallchains_hom1chain', '*.ssv')))

    # rp_raw_cifs = sorted(glob.glob(os.path.join(_rp_nmr_dir(), 'raw_cifs', 'hethom_combined', '*.cif')))
    # for i, bla in enumerate(rp_raw_cifs):
    #     if os.path.basename(bla).removesuffix('.cif') == '8ONU':
    #         print(i)
    # pass
    # ca_counts_ = _calc_ca_counts(rp_parsed_cifs_ssvs_)

    # ca_counts_dir = os.path.join(_rp_nmr_dir(), 'stats', 'multimod_2713_hetallchains_hom1chain')
    # os.makedirs(ca_counts_dir, exist_ok=True)
    # ca_counts_lst_f = os.path.join(ca_counts_dir, 'ca_counts.lst')

    # with open(ca_counts_lst_f, 'w') as f:
    #     f.writelines(str(s) + '\n' for s in ca_counts_)
    # plot_protein_lengths(ca_counts_)
    # plot_protein_lengths_binned(ca_counts_)
    # with open(ca_counts_lst_f, 'r') as f:
    #     ca_counts_ = f.read().splitlines()
    # ca_counts_ = [int(ca_count_) for ca_count_ in ca_counts_]
    # ca_counts_.sort()
    # plot_protein_lengths_binned(ca_counts_, bin_size=10)

    rp_pidChains_lst_f_ = os.path.join('..', 'data', 'NMR', 'multimodel_lists', 'multimod_2713_hetallchains_hom1chain.lst')

    stats_pdf = generate_stats(sub_dir='multimod_2713_hetallchains_hom1chain',
                               rp_pidchains_lst_f=rp_pidChains_lst_f_,
                               rp_fasta_f= rp_fasta_f_,
                               run_and_write_mmseqs2=False,
                               run_and_write_rmsd=False, use_mmcif=True)
    stats_dst_dir = os.path.join(_rp_nmr_dir(), 'stats', 'multimod_2713_hetallchains_hom1chain')
    os.makedirs(stats_dst_dir, exist_ok=True)
    stats_dst_f = os.path.join(stats_dst_dir, 'multimod_2713_hetallchains.csv')
    stats_pdf.to_csv(stats_dst_f, index=False)

    # TODO generate 2d plots (akin to bar graphs) of datasets for:
    #   - protein lengths (CA count)
    #   - rmsd min and max values
    #   - number of models
    #   - date of deposition

