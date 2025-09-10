import os, glob, math, json
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import RMSD

# RELATIVE PATHS:
def _rp_nmr_dir() -> str:
    return os.path.join('..', 'data', 'NMR')

def _rp_matrices_dir(sub_dir: str, rmsd_or_tms: str) -> str:
    rmsd_or_tms.lower()
    dir_ = 'RMSD' if rmsd_or_tms == 'rmsd' else 'TMscores'
    return os.path.join(_rp_nmr_dir(), dir_, sub_dir, f'{rmsd_or_tms}_matrices')

def _rp_stats_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'stats', sub_dir)

def plot_counts(attribute: str, attr_counts: list, bin_size: int=1, num_pidchains=2713):
    attr_counts = np.array(attr_counts)

    if bin_size == 1:
        heights = attr_counts
        x = np.arange(len(attr_counts))
    else:
        n_bins = len(attr_counts) // bin_size
        trimmed = attr_counts[:n_bins * bin_size]
        grouped = trimmed.reshape(n_bins, bin_size)
        heights = grouped.mean(axis=1)
        x = np.arange(n_bins)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x, heights, color='lightgrey', edgecolor='gainsboro', linewidth=0.5, width=1.0, align='edge')
    ax.set_xlim(left=-5)

    ax.set_xlabel(f'PDBchains (bin size={bin_size})', fontsize=10)

    target_n_ticks = 30
    tick_interval = max(1, int(math.ceil(len(x) / target_n_ticks / 10.0)) * 10)
    tick_positions = np.arange(0, len(x) + 1, tick_interval)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_positions, rotation=-90, fontsize=8)

    ax.set_ylabel(f'{attribute} count', fontsize=10)
    ax.set_title(f'Number of {attribute} in each PDBchain, (for {num_pidchains} PDBchains)', fontsize=12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('lightgrey')
    ax.spines['bottom'].set_color('lightgrey')
    fig.tight_layout()
    plt.show()

def _calc_model_counts(rp_parsed_cifs_ssvs: list) -> list:
    model_counts = list()
    pidc_list, model_count_list = [], []
    model_counts_pidc = {'pidc': pidc_list, 'model_counts': model_count_list}
    for rp_parsed_cif_ssv in rp_parsed_cifs_ssvs:
        pdf = pd.read_csv(rp_parsed_cif_ssv, sep=' ')
        model_count = len(pdf['A_pdbx_PDB_model_num'].unique())
        pidc = os.path.basename(rp_parsed_cif_ssv).removesuffix('.ssv')
        pidc_list.append(pidc)
        model_count_list.append(model_count)
        model_counts.append(model_count)
    model_counts_pidc['pidc'] = pidc_list
    model_counts_pidc['model_counts'] = model_count_list
    model_counts_pidc_pdf = pd.DataFrame(model_counts_pidc)
    model_counts_pidc_pdf = model_counts_pidc_pdf.sort_values(by=['model_counts'], ascending=[True])
    rp_dst_csv = os.path.join(_rp_stats_dir('multimod_2713_hetallchains_hom1chain'), 'model_counts.csv')
    model_counts_pidc_pdf.to_csv(rp_dst_csv, index=False)
    return sorted(model_counts)

def _calc_ca_counts(rp_parsed_cifs_ssvs: list) -> list:
    ca_counts = list()
    pidc_list, ca_count_list = [], []
    ca_counts_pidc = {'pidc': pidc_list, 'ca_counts': ca_count_list}

    for rp_parsed_cif_ssv in rp_parsed_cifs_ssvs:
        pdf = pd.read_csv(rp_parsed_cif_ssv, sep=' ')
        ca_counts_all_models = pdf.shape[0]
        model_count = len(pdf['A_pdbx_PDB_model_num'].unique())
        ca_count = int(ca_counts_all_models / model_count)
        ca_counts.append(ca_count)
        pidc = os.path.basename(rp_parsed_cif_ssv).removesuffix('.ssv')
        pidc_list.append(pidc)
        ca_count_list.append(ca_count)
        if ca_count < 3:
            print(f'{pidc} has < 3 CAs, with only {ca_count} CAs.')
    ca_counts_pidc['pidc'] = pidc_list
    ca_counts_pidc['ca_counts'] = ca_count_list
    ca_counts_pidc_pdf = pd.DataFrame(ca_counts_pidc)
    ca_counts_pidc_pdf = ca_counts_pidc_pdf.sort_values(by=['ca_counts'], ascending=[True])
    rp_dst_csv = os.path.join(_rp_stats_dir('multimod_2713_hetallchains_hom1chain'), 'ca_counts.csv')
    ca_counts_pidc_pdf.to_csv(rp_dst_csv, index=False)
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

def violin_plot(pdf):
    plt.figure(figsize=(14, 6))
    sns.violinplot(data=pdf, x='year', y='total_model_count', inner='quartile', color='skyblue', linewidth=1)
    sns.swarmplot(data=pdf, x='year', y='total_model_count', size=3, color='black', alpha=0.6)
    plt.title('Number of model by deposition year', fontsize=14)
    plt.xlabel('Year', fontsize=12)
    plt.ylabel('Model count', fontsize=12)
    plt.xticks(rotation=-90)
    plt.grid(True, linestyle='--', alpha=0.3)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.show()

def plot_rmsds_and_stdev(pdf):
    # Sort by mean_rmsd (optional: for more readable axis)
    pdf = pdf.sort_values('mean_rmsd').reset_index(drop=True)

    x = range(len(pdf))
    y = pdf['mean_rmsd']
    yerr = pdf['stdev_rmsd']
    ymin = pdf['min_rmsd']
    ymax = pdf['max_rmsd']
    labels = pdf['Pid_chain']

    fig, ax = plt.subplots(figsize=(14, 6))

    # Draw vertical lines from min to max rmsd
    ax.vlines(x, ymin, ymax, color='lightgrey', alpha=0.7, linewidth=1)

    # Overlay mean RMSD points with SEM error bars
    ax.errorbar(x, y, yerr=yerr, fmt='o', color='steelblue', markersize=3, capsize=2, linewidth=1)

    # Clean up axis
    ax.set_xlabel('PDB Chain', fontsize=10)
    ax.set_ylabel('RMSD', fontsize=10)
    ax.set_title('RMSD per PDB Chain (mean ± SEM, min/max range)', fontsize=12)
    ax.set_xlim(-1, len(pdf))  # pad edges
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, linestyle='--', alpha=0.3)

    # Optional: thin x-ticks for readability
    ax.set_xticks(x[::200])  # show every 200th label only
    ax.set_xticklabels(labels[::200], rotation=-90, fontsize=6)

    plt.tight_layout()
    plt.show()

def plot_rmsd_and_tms(pdf):
    # sns.set_theme(style='whitegrid', context='talk')
    fig, ax1 = plt.subplots(figsize=(8, 5))

    pidc = pdf['Pidc'].values
    mean_tms_vals = pdf['meanTMS'].values
    mean_rmsd_vals = pdf['meanRMSD'].values

    ax1.set_xlabel('PDB-chain')
    ax1.set_ylabel('mean TM-score', color='tab:blue', fontsize=10)
    sns.scatterplot(x=pidc, y=mean_tms_vals, ax=ax1, color='tab:blue', s=10)
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.spines['top'].set_visible(False)
    step = 50
    ax1.set_xticks(np.arange(0, len(pidc), step))
    # ax1.set_xticklabels(np.arange(0,  + 1, 20), rotation=-90, fontsize=8)
    ax1.set_xticklabels(pidc[::step], rotation=-90, fontsize=6)


    ax2 = ax1.twinx()
    ax2.spines['top'].set_visible(False)
    ax2.set_ylabel('mean RMSD', color='tab:red', fontsize=10, rotation=-90, labelpad=15)
    sns.scatterplot(x=pidc, y=mean_rmsd_vals, ax=ax2, color='tab:red', s=10)
    ax2.tick_params(axis='y', labelcolor='tab:red')


    plt.title('Compare mean TM-scores to mean RMSDs for each PDB-chain')
    plt.tight_layout()
    plt.show()

def dendrgm_heatmap_contourmap(dist_matrix, pidc: str, rmsd_or_tms: str):
    n_models = dist_matrix.shape[0]
    linkage_mat = RMSD.dendrogrm(dist_matrix, n_models, pidc, rmsd_or_tms)
    RMSD.heatmap(dist_matrix, linkage_mat)
    # RMSD.contour_map(dist_matrix)


if __name__ == '__main__':
    rmsd_or_tms_ = 'tms'
    # rmsd_or_tms_ = 'rmsd'
    sub_dir_ = 'multimod_2713_hetallchains_hom1chain'
    rp_matrices_dir = os.path.join(_rp_matrices_dir(sub_dir_, rmsd_or_tms_))
    rp_matrices = sorted(glob.glob(os.path.join(rp_matrices_dir, '*.npz')))
    results, res = list(), dict()
    for rp_matrix in rp_matrices[1:50]:
        pidc = os.path.basename(rp_matrix).removesuffix('.npz')
        print(pidc)
        loaded = np.load(rp_matrix)
        dist_matrix = loaded['mat']
        if rmsd_or_tms_.lower() == 'tms':  # need to convert similarity matrix to distance matrix
            dist_matrix = 1.0 - dist_matrix
        dendrgm_heatmap_contourmap(dist_matrix, pidc, rmsd_or_tms_)
        if rmsd_or_tms_ == 'tms':
            thrshld = 0.5
        else:
            thrshld = 6
        clusters_dict = RMSD.cluster_models(dist_matrix, threshold=thrshld)   # WHAT THRESHOLD TO USE ??
        rp_clusters_dir = os.path.join(rp_matrices_dir, f'clusters_{str(thrshld).replace('.', 'p')}')
        os.makedirs(rp_clusters_dir, exist_ok=True)

        clusters_dict = {int(k): v for k, v in clusters_dict.items()}

        with open(os.path.join(rp_clusters_dir, f'{pidc}.json'), 'w') as f:
            json.dump(clusters_dict, f)


    # rp_rmsd_mat_dir = os.path.join(_rp_rmsd_dir('multimod_2713_hetallchains_hom1chain'), 'rmsd_matrices')
    # rp_rmsd_mats = sorted(glob.glob(os.path.join(rp_rmsd_mat_dir, '*.npz')))
    # for rp_rmsd_mat in rp_rmsd_mats[:10]:
    #     dendrgm_heatmap_contourmap(rp_rmsd_mat)
    # pass
    # rp_stats_pidc_2713_csv = os.path.join('..', 'data', 'NMR', 'stats',
    #                                  'multimod_2713_hetallchains_hom1chain', 'mm_2713_pidc.csv')
    # pidc_2713_pdf = pd.read_csv(rp_stats_pidc_2713_csv)
    # pidc_2713_pdf = pidc_2713_pdf.sort_values(by=['meanRMSD'], ascending=False)
    # plot_rmsd_and_tms(pidc_2713_pdf)
