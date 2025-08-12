import os, glob
import statistics as stats
from time import time
from typing import Tuple
import math
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio import SeqIO
import RMSD
import mmseqs2
import tm_aligner
from src.tm_aligner import rp_tmscores_dir


# BUILDING RELATIVE PATHS:
def _rp_nmr_dir() -> str:
    return os.path.join('..', 'data', 'NMR')

def _rp_stats_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'stats', sub_dir)

def _rp_rmsd_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'RMSD', sub_dir)

def _rp_pdb_chains_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'pdb_chains', sub_dir)

def _rp_parsed_cifs_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', sub_dir)

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


def pdbid_dict_chain(pdbid_chains: list) -> dict:
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
    return dict(pdbid_chains_dict)


def _read_multimodel_pdbid_chains(rp_pidChains: str) -> Tuple[list, dict, list]:

    if rp_pidChains.endswith('.txt'):
        with open(rp_pidChains, 'r') as f:
            # pdbid_chains = [line.rstrip('\n') for line in f]
            pdbid_chains = f.read().split()
    else:  # .lst file
        with open(rp_pidChains, 'r') as f:
            pdbid_chains = f.read().splitlines()

    pdbid_chains.sort()
    pidchains_dict = pdbid_dict_chain(pdbid_chains)
    pdbids = list(pidchains_dict.keys())
    return pdbids, pidchains_dict, pdbid_chains


def _calc_identity_for_stats(het_hom: str, pdbid: str, filt_pdf) -> str:
    if not filt_pdf.empty:
        print(f'There is some heterogeneity of sequence between some chains. '
              f'See data/NMR/mmseqs/{het_hom}/results/{pdbid}.csv')
        identity = '(chains not identical, see mmseqs pdf)'
    else:
        print('All chains have high identity, if not 100%')
        identity = '0'
    return identity


def _compute_rmsd(run_and_write_rmsd: bool, pidc: str, sub_dir: str):
    print(pidc)
    rp_parsed_cif_dir = _rp_parsed_cifs_dir(sub_dir)
    rp_parsed_cifs_ssv = os.path.join(rp_parsed_cif_dir, f'{pidc}.ssv')
    rp_rmsd_mean_coords_dir = os.path.join(_rp_rmsd_dir(sub_dir), 'mean_coords')
    rp_mean_coords_csv = os.path.join(rp_rmsd_mean_coords_dir, f'{pidc}.csv')
    # GENERATE NEW RMSD SCORES OR READ PRE-WRITTEN VALUES:
    if run_and_write_rmsd:
        rmsds, model_nums, pdc4pdf = RMSD.calc_rmsds_of_models(rp_parsed_cifs_ssv, rp_mean_coords_csv)
        rmsd_pdf = pd.DataFrame({'pidc': pdc4pdf, 'model_num': model_nums, 'rmsd': rmsds})
        rp_rmsd_dir = os.path.join(_rp_rmsd_dir(sub_dir))
        os.makedirs(rp_rmsd_dir, exist_ok=True)
        rp_rmsd_csv = os.path.join(rp_rmsd_dir, f'{pidc}.csv')
        rmsd_pdf.to_csv(rp_rmsd_csv, index=False)
    else:  # This relies on it having been run separately beforehand, and we can just read it in now:
        rp_rmsd_per_model_dir = os.path.join(_rp_rmsd_dir(sub_dir))
        rp_rmsd_per_model_csv_f = os.path.join(rp_rmsd_per_model_dir, f'{pidc}.csv')
        rmsd_pdf = pd.read_csv(rp_rmsd_per_model_csv_f)
    rmsds = rmsd_pdf['rmsd'].values
    if len(rmsds) > 0:
        min_rmsd, max_rmsd, mean_rmsd, stdev_rmsd = np.min(rmsds), np.max(rmsds), stats.mean(rmsds), stats.stdev(rmsds)
        # sem_rmsd = stdev_rmsd / np.sqrt(len(rmsds))
    else:
        min_rmsd, max_rmsd, mean_rmsd, stdev_rmsd = np.nan, np.nan, np.nan, np.nan
    return min_rmsd, max_rmsd, mean_rmsd, stdev_rmsd

def _compute_tms(run_and_write_tms: bool, pidc: str, sub_dir: str) -> tuple:
    
    if not run_and_write_tms: # THEN READ IN ALREADY CALCULATED VALUES
        rp_tmscores_dir = tm_aligner.rp_tmscores_dir(sub_dir)
        rp_tms_csv_f = os.path.join(rp_tmscores_dir, f'TMS_{pidc}.csv')
        pidc_tms_pdf = pd.read_csv(rp_tms_csv_f)
        min_tms = pidc_tms_pdf['min_TMS'].values[0]
        max_tms = pidc_tms_pdf['max_TMS'].values[0]
        mean_tms = pidc_tms_pdf['mean_TMS'].values[0]
        stdev_tms = pidc_tms_pdf['stdev_TMS'].values[0]
    else:
        # NOTE: TM-ALIGN RETURNS NONE FOR PROTEINS WITH LESS THAN 3 CA ATOMS (WHICH THESE 11 PDB-CHAINS HAVE):
        if (pidc == '1GAC_A' or pidc == '1GAC_B' or pidc == '1GAC_C' or pidc == '1GAC_D' or
                pidc == '1WCO_A' or pidc == '2AIZ_B' or pidc == '2K1Q_B' or pidc == '2M9P_B' or
                pidc == '2M9Q_B' or pidc == '2MX6_B' or pidc == '3CYS_B'):
            min_tms, max_tms, mean_tms, stdev_tms = np.nan, np.nan, np.nan, np.nan

        rp_pidc_2713_dir = _rp_pdb_chains_dir('multimod_2713_hetallchains_hom1chain')
        rp_mean_coords_dir = os.path.join(rp_pidc_2713_dir, 'mean_coords')
        rp_mean_coords_pdb = os.path.join(rp_mean_coords_dir, f'{pidc}.pdb')
        rp_permodel_pdb_pidc_dir = os.path.join(rp_pidc_2713_dir, 'per_model', pidc)
        tm_scores = []
        rp_permodel_pdbs = sorted(glob.glob(os.path.join(rp_permodel_pdb_pidc_dir, '*.pdb')))
        for rp_permodel_pdb in rp_permodel_pdbs:
            tms = tm_aligner.compute_tm(rp_pdb1=rp_permodel_pdb, rp_pdb2=rp_mean_coords_pdb)
            tm_scores.append(tms)
        min_tms, max_tms, mean_tms, stdev_tms = min(tm_scores), max(tm_scores), np.mean(tm_scores), np.std(tm_scores)
    return min_tms, max_tms, mean_tms, stdev_tms



def generate_stats(sub_dir: str, rp_pidc_lst_f: str, rp_fasta_f: str, run_and_write_mmseqs2=False,
                   run_and_write_rmsd=False, run_and_write_tms=False, use_mmcif=True) -> tuple:
    """
    (Takes 18 mins to complete 2713 PDB-chains.)
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
    pid_stats = {'Pid': [], 'chains': [], 'hethom': [], 'year': [], 'protchains_counts': []}
    pids_list, chains_list, hethoms_list, years_list, protchains_counts_list= [], [], [], [], []

    pidc_stats = {'Pidc': [], 'CA_count': [],
                  'minRMSD': [], 'maxRMSD': [], 'meanRMSD': [], 'stdevRMSD': [],
                  'minTMS': [], 'maxTMS': [], 'meanTMS': [], 'stdevTMS': [],
                  'year': [], 'models_count': [], 'homologues': []}
    pidc_list, ca_counts_list = [], []
    min_rmsd_list, max_rmsd_list, mean_rmsd_list, stdev_rmsd_list = [], [], [], []
    min_tms_list, max_tms_list, mean_tms_list, stdev_tms_list = [], [], [], []
    years_list, models_counts_list, homologues_list = [], [], []

    fname = os.path.basename(rp_fasta_f).removesuffix('.fasta')
    rp_mmseqs2_results_dir = os.path.join(_rp_mmseqs_dir(sub_dir), 'results')
    os.makedirs(rp_mmseqs2_results_dir, exist_ok=True)
    rp_mmseqs2_csv_f = os.path.join(rp_mmseqs2_results_dir, f'{fname}.csv')

    if use_mmcif:
        pdb_or_cif = 'cif'
        parser = MMCIFParser(QUIET=True)  # Used in _total_chain_count_and_year() below (~line 150).
        # 'QUIET' is for debugging any problems that might pop up in PDB records.
    else:
        pdb_or_cif = 'pdb'
        parser = PDBParser(QUIET=True)

    # NOTE: Next 7 lines ARE ONLY FOR assigning het/hom (~line 199) in stats pdf,
    # and for assertion that each raw file available (~line 159):
    rp_raw_struct_het_dir = os.path.join(_rp_nmr_dir(), f'raw_{pdb_or_cif}s', 'heteromeric')
    rp_raw_struct_hom_dir = os.path.join(_rp_nmr_dir(), f'raw_{pdb_or_cif}s', 'homomeric')
    rp_raw_het_struct_files = glob.glob(os.path.join(rp_raw_struct_het_dir, f'*.{pdb_or_cif}'))
    rp_raw_hom_struct_files = glob.glob(os.path.join(rp_raw_struct_hom_dir, f'*.{pdb_or_cif}'))
    rp_raw_struct_files = rp_raw_het_struct_files + rp_raw_hom_struct_files
    raw_pids = [os.path.basename(rp_raw_struct_f).removesuffix(f'.{pdb_or_cif}') for rp_raw_struct_f in rp_raw_struct_files]
    raw_het_pids = [os.path.basename(rp_raw_het_struct_f).removesuffix(f'.{pdb_or_cif}')
                    for rp_raw_het_struct_f in rp_raw_het_struct_files]

    # GENERATE NEW MMSEQS2 ALIGNMENT SCORES OR READ PRE-WRITTEN VALUES:
    if run_and_write_mmseqs2:
        mmseq2_pdf = mmseqs2.run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f)
        mmseq2_pdf = mmseqs2.dedupe_rm_self_aligns(mmseq2_pdf)
        mmseq2_pdf.to_csv(rp_mmseqs2_csv_f, index=False)  # shape=(33074, 7)
    else:
        mmseq2_pdf = pd.read_csv(rp_mmseqs2_csv_f)  # shape=(33074, 7)

    pids, pid_chains_dict, pidc_list = _read_multimodel_pdbid_chains(rp_pidc_lst_f)
    rp_raw_struct_dir = os.path.join(_rp_nmr_dir(), f'raw_{pdb_or_cif}s', 'hethom_combined')
    rp_parsed_cifs_dir = os.path.join(_rp_parsed_cifs_dir(sub_dir))

    for i, (pid, chains) in enumerate(pid_chains_dict.items()):
        if use_mmcif:
            assert pid in raw_pids, f'{pid}.cif not round in raw_cifs.'
        else:
            if pid not in ['9D9B', '7ZE0', '9D9C', '9D9A']:  # (These 4 are not found in 'legacy' PDB files on RCSB)
                assert pid in raw_pids, f'{pid}.pdb not round in raw_pdbs.'
            else:
                print(f'{pid}.pdb is not available in legacy PDB for downloading from RCSB. So, cannot include it.')
                continue
        rp_raw_struct_f = os.path.join(rp_raw_struct_dir, f'{pid}.{pdb_or_cif}')
        total_chain_count, year = _total_chain_count_and_year(rp_raw_struct_f, parser)
        het_hom = 'het' if pid in raw_het_pids else 'hom'
        pids_list.append(pid)
        chains_list.append(chains)
        hethoms_list.append(het_hom)
        years_list.append(year)
        protchains_counts_list.append(total_chain_count)

        for chain in chains:
            pidc = f'{pid}_{chain}'
            pidc_list.append(pidc)

            min_rmsd, max_rmsd, mean_rmsd, stdev_rmsd = _compute_rmsd(run_and_write_rmsd, pidc, sub_dir)
            min_rmsd_list.append(min_rmsd)
            max_rmsd_list.append(max_rmsd)
            mean_rmsd_list.append(mean_rmsd)
            stdev_rmsd_list.append(stdev_rmsd)

            min_tms, max_tms, mean_tms, stdev_tms = _compute_tms(run_and_write_tms, pidc, sub_dir)
            min_tms_list.append(min_tms)
            max_tms_list.append(max_tms)
            mean_tms_list.append(mean_tms)
            stdev_tms_list.append(stdev_tms)

            # GET MODELS_COUNT, CA_COUNT FROM PARSED CIF SSV:
            rp_parsed_cif_ssv_f = os.path.join(rp_parsed_cifs_dir, f'{pidc}.ssv')
            pidc_pdf = pd.read_csv(rp_parsed_cif_ssv_f, sep=' ')
            models_count = len(pidc_pdf['A_pdbx_PDB_model_num'].unique())
            models_counts_list.append(models_count)
            ca_counts_list.append(int(pidc_pdf.shape[0] / models_count))

            # GET LIST OF HOMOLOGOUS PDB-CHAHINS FOR GIVEN PDB-CHAIN:
            homologues = mmseqs2.find_homologues_30_20_90(mmseq2_pdf, pidc)
            homologues_list.append(homologues)

        pidc_stats['Pidc'] = pidc_list,
        pidc_stats['CA_count'] = ca_counts_list,
        pidc_stats['minRMSD'] = min_rmsd_list,
        pidc_stats['maxRMSD'] = max_rmsd_list,
        pidc_stats['meanRMSD'] = mean_rmsd_list,
        pidc_stats['stdevRMSD'] = stdev_rmsd_list,
        pidc_stats['minTMS'] = min_tms_list,
        pidc_stats['maxTMS'] = max_tms_list,
        pidc_stats['meanTMS'] = mean_tms_list,
        pidc_stats['stdevTMS'] = stdev_tms_list,
        pidc_stats['year'] = years_list,
        pidc_stats['models_count'] = models_counts_list,
        pidc_stats['homologues'] = homologues_list

    pid_stats['Pid'] = pids_list
    pid_stats['chains'] = chains_list
    pid_stats['hethom'] = hethoms_list
    pid_stats['year'] = years_list
    pid_stats['protchains_counts'] = protchains_counts_list
    pid_pdf = pd.DataFrame(pid_stats)
    pid_pdf = pid_pdf.sort_values(by=['year'], ascending=[True])

    pidc_pdf = pd.DataFrame(pidc_stats)
    pidc_pdf = pidc_pdf.sort_values(by=['maxRMSD', 'maxTMS'], ascending=[False, True])

    print(pidc_pdf.dtypes)
    for col_to_cast in ['Pid', 'chain']:
        pidc_pdf[col_to_cast] = pidc_pdf[col_to_cast].astype('string')
    print(pidc_pdf.dtypes)

    print(f'Completed {len(pidc_list)} PDB-chains in {round((time() - start) / 60)} mins')
    return pid_pdf, pidc_pdf


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


# COPIED OVER TO RMSD.py
def _calc_rmsds_stats(pidchains: list):
    """
    # NOTE: COPIED OVER TO RMSD.py
    TODO Note that rmsds for 6UJV_A, 7CLV_A, 7CLV_B & 8J4I_A are empty.. need to have a closer look at this to see why...
    """
    rmsd_stats_per_pidc = []
    for pid_chain in pidchains:
        rp_rmsd_per_model_dir = os.path.join(_rp_nmr_dir(), 'RMSD', 'multimod_2713_hetallchains_hom1chain')
        rp_rmsd_per_model_csv_f = os.path.join(rp_rmsd_per_model_dir, f'{pid_chain}.csv')
        rmsd_pdf = pd.read_csv(rp_rmsd_per_model_csv_f)
        rmsds = rmsd_pdf['rmsd'].values
        rp_parsed_cifs_dir = os.path.join(_rp_parsed_cifs_dir(), 'multimod_2713_hetallchains_hom1chain')
        pdf = pd.read_csv(os.path.join(rp_parsed_cifs_dir, f'{pid_chain}.ssv'), sep=' ')
        model_count = len(pdf['A_pdbx_PDB_model_num'].unique())

        if len(rmsds) > 0:
            min_rmsd, max_rmsd = np.min(rmsds), np.max(rmsds)
            mean_rmsd, stdev_rmsd = np.mean(rmsds), np.std(rmsds)
        else:
            print(f'rmsds for {pid_chain} is empty: {rmsds}. Cannot calc min/max. Cannot include.')
            continue
        rmsd_stats_per_pidc.append({
            'Pid_chain': pid_chain,
            'min_rmsd': min_rmsd,
            'max_rmsd': max_rmsd,
            'mean_rmsd': mean_rmsd,
            'stdev_rmsd': stdev_rmsd,
            'model_count': model_count,
        })
    rmsdstats_pdf = pd.DataFrame(rmsd_stats_per_pidc)
    rmsdstats_pdf = rmsdstats_pdf.sort_values(by=['mean_rmsd', 'stdev_rmsd'], ascending=[True, True])
    return rmsdstats_pdf


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


def _tabulate_year_chain_model_counts(rp_raw_cif_files: list):
    start = time()
    rp_raw_cifs_het_dir = os.path.join(_rp_nmr_dir(), f'raw_cifs', 'heteromeric')
    rp_raw_het_cifs_files = glob.glob(os.path.join(rp_raw_cifs_het_dir, f'*.cif'))
    het_pids = [os.path.basename(rp_raw_het_cif_f).removesuffix(f'.cif')
                    for rp_raw_het_cif_f in rp_raw_het_cifs_files]
    biopy_parser = MMCIFParser(QUIET=True)
    year_chain_model_counts = []
    for rp_raw_cif_f in rp_raw_cif_files:
        bio_struct = biopy_parser.get_structure('', rp_raw_cif_f)
        total_model_count = len(list(bio_struct.get_models()))
        total_chain_count = len(list(next(bio_struct.get_models()).get_chains()))
        year = bio_struct.header['deposition_date'][:4]
        pid = os.path.basename(rp_raw_cif_f).removesuffix('.cif')
        print(pid)
        het_hom = 'het' if pid in het_pids else 'hom'

        year_chain_model_counts.append({
            'PDBid': pid,
            'het_hom': het_hom,
            'year': int(year),
            'total_model_count': int(total_model_count),
            'total_chain_count': int(total_chain_count)
        })
    ycmc_pdf = pd.DataFrame(year_chain_model_counts)
    ycmc_pdf = ycmc_pdf.sort_values(by=['year', 'total_model_count'], ascending=[True, True])
    print(f'Completed {len(rp_raw_cif_files)} PDBs in {round((time() - start) / 60)} mins')
    return ycmc_pdf


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
    for rp_parsed_cif_ssv in rp_parsed_cifs_ssvs:
        pdf = pd.read_csv(rp_parsed_cif_ssv, sep=' ')
        model_count = len(pdf['A_pdbx_PDB_model_num'].unique())
        if model_count == 1:
            print(f'Single model PDB: {os.path.basename(rp_parsed_cif_ssv).removesuffix('.ssv')}')
        model_counts.append(model_count)
    return sorted(model_counts)


def _calc_ca_counts(rp_parsed_cifs_ssvs: list) -> list:
    ca_counts = list()
    ca_counts_pidc = {'pidc': [], 'ca_counts': []}
    for rp_parsed_cif_ssv in rp_parsed_cifs_ssvs:
        pdf = pd.read_csv(rp_parsed_cif_ssv, sep=' ')
        ca_counts_all_models = pdf.shape[0]
        model_count = len(pdf['A_pdbx_PDB_model_num'].unique())
        ca_count = int(ca_counts_all_models / model_count)
        ca_counts.append(ca_count)

        pidc = os.path.basename(rp_parsed_cif_ssv).removesuffix('.ssv')
        ca_counts_pidc['pidc'].append(pidc)
        ca_counts_pidc['ca_counts'].append(ca_count)
        if ca_count <= 3:
            print(f'{pidc} has <4 CAs, with only {ca_count} CAs.')

    ca_counts_pidc_pdf = pd.DataFrame(ca_counts_pidc)
    ca_counts_pidc_pdf = ca_counts_pidc_pdf.sort_values(by=['ca_counts'], ascending=[True])
    rp_dst_csv = os.path.join(_rp_stats_dir('multimod_2713_hetallchains_hom1chain'), 'model_counts.csv')
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


if __name__ == '__main__':

    # rp_raw_cifs = sorted(glob.glob(os.path.join(_rp_nmr_dir(), 'raw_cifs', 'hethom_combined', '*.cif')))

    # # VISUALISATIONS:

    # # 1. FASTA SEQUENCE DISTRIBUTION:
    # rp_fasta_f_ = os.path.join(_rp_mmseqs_fasta_dir(sub_dir='multimod_2713_hetallchains_hom1chain'),
    #                            'multimod_2713_hetallchains_hom1chain.fasta')
    # plot_fasta_size_distribution(rp_fasta_f_, x_limit_220=False)

    # # 2. BAR CHART OF ALPHA-CARBON COUNT FOR EACH PDB:
    # # 2.A. CALCULATE CA COUNT FROM PARSED CIFS DATA AND WRITE TO STATS/...LST FILE:
    rp_parsed_cifs_ssvs_ = sorted(glob.glob(os.path.join(_rp_parsed_cifs_dir('multimod_2713_hetallchains_hom1chain'),
                                                         '*.ssv')))
    ca_counts_ = _calc_ca_counts(rp_parsed_cifs_ssvs_)
    # ca_counts_dir = _rp_stats_dir('multimod_2713_hetallchains_hom1chain')
    # os.makedirs(ca_counts_dir, exist_ok=True)
    # ca_counts_lst_f = os.path.join(ca_counts_dir, 'ca_counts.lst')
    # with open(ca_counts_lst_f, 'w') as f:
    #     f.writelines(str(s) + '\n' for s in ca_counts_)

    # # 2.B. READ CA COUNT FROM STATS/...LST FILE AND PLOT CA COUNT:
    # ca_counts_dir = _rp_stats_dir('multimod_2713_hetallchains_hom1chain')
    # ca_counts_lst_f = os.path.join(ca_counts_dir, 'ca_counts.lst')
    # with open(ca_counts_lst_f, 'r') as f:
    #     ca_counts_ = f.read().splitlines()
    # ca_counts_ = [int(ca_count_) for ca_count_ in ca_counts_]
    # ca_counts_.sort()
    # plot_counts(attribute='CA', attr_counts=ca_counts_, bin_size=10, num_pidchains=2713)

    # # 3. BAR CHART OF MODEL COUNT FOR EACH PDB:
    # # 3.A. CALCULATE MODEL COUNT FROM PARSED CIFS DATA AND WRITE TO STATS/...LST FILE:
    # rp_parsed_cifs_ssvs_ = sorted(glob.glob(os.path.join(_rp_parsed_cifs_dir(),
    #                                              'multimod_2713_hetallchains_hom1chain', '*.ssv')))
    # model_counts_ = _calc_model_counts(rp_parsed_cifs_ssvs_)
    # model_counts_dir = _rp_stats_dir('multimod_2713_hetallchains_hom1chain')
    # os.makedirs(model_counts_dir, exist_ok=True)
    # model_counts_lst_f = os.path.join(model_counts_dir, 'model_counts.lst')
    # with open(model_counts_lst_f, 'w') as f:
    #     f.writelines(str(s) + '\n' for s in model_counts_)

    # # 3.B. READ MODEL COUNT FROM STATS/...LST FILE AND PLOT MODEL COUNT:
    # model_counts_dir = _rp_stats_dir('multimod_2713_hetallchains_hom1chain'
    # model_counts_lst_f = os.path.join(model_counts_dir, 'model_counts.lst')
    # with open(model_counts_lst_f, 'r') as f:
    #     model_counts_ = f.read().splitlines()
    # model_counts_ = [int(model_count_) for model_count_ in model_counts_]
    # model_counts_.sort()
    # plot_counts(attribute='model', attr_counts=model_counts_, bin_size=10, num_pidchains=2713)

    # # 4. YEAR (DEPOSITION), MODEL COUNT, CHAIN COUNT - FROM ALL DOWNLOADED NMR CIF FILES (NOT THE PDBCHAINS):
    # # 4.A. EXTRACT YEAR, MODEL COUNT, CHAIN COUNT, FROM PDB FILE AND WRITE TO STATS/...CSV:
    # rp_raw_cif_files = sorted(glob.glob(os.path.join(_rp_nmr_dir(), 'raw_cifs', 'hethom_combined', '*.cif')))
    # ycmc_pdf_ = _tabulate_year_chain_model_counts(rp_raw_cif_files)
    # year_modcounts_dir = _rp_stats_dir('from_1725_raw_cifs')
    # os.makedirs(year_modcounts_dir, exist_ok=True)
    # year_modcounts_csv_f = os.path.join(year_modcounts_dir, 'year_chain_model_counts.csv')
    # ycmc_pdf_.to_csv(year_modcounts_csv_f, index=False)

    # # 4.B. READ CSV FOR EXTRACT YEAR, MODEL COUNT, CHAIN COUNT AND PLOT:
    # year_modcounts_dir = _rp_stats_dir('from_1725_raw_cifs')
    # year_modcounts_csv_f = os.path.join(year_modcounts_dir, 'year_chain_model_counts.csv')
    # ycmc_pdf_ = pd.read_csv(year_modcounts_csv_f)
    # violin_plot(ycmc_pdf_)

    # # 5. MIN, MAX, MEAN & STDDEV OF RMSD VALUES FOR EACH PDBCHAIN'S MODEL AGAINST THE MEAN COORDS OF ALL MODELS.
    # # 5.A. CALCULATE VALUES AND WRITE TO STATS/...CSV:
    # rp_pidchains_lst_f_ = os.path.join(_rp_nmr_dir(), 'multimodel_lists',
    #                                    'multimod_2713_hetallchains_hom1chain.lst')
    # with open(rp_pidchains_lst_f_, 'r') as f:
    #     pidchains_2713 = sorted(f.read().splitlines())
    # rsmds_stats_pdf = _calc_rmsds_stats(pidchains_2713)
    # rp_rsmds_stats_csv_f = os.path.join(_rp_stats_dir('multimod_2713_hetallchains_hom1chain'), 'rmsds_stats.csv')
    # rsmds_stats_pdf.to_csv(rp_rsmds_stats_csv_f, index=False)

    # # 5.B. READ CSV FOR MIN, MAX, MEAN & STDDEV RSMD VALUES FOR 2713 PDB-CHAINS:
    # rp_rsmds_stats_csv_f = os.path.join(_rp_stats_dir('multimod_2713_hetallchains_hom1chain'), 'rmsds_stats.csv')
    # rsmds_stats_pdf = pd.read_csv(rp_rsmds_stats_csv_f)
    # plot_rmsds_and_stdev(rsmds_stats_pdf)
    # pass

    # # MAIN STATS FUNCTION
    # # GENERATE STATS PDF AND WRITE TO CSV: (Takes 18 mins to complete 2713 PDB-chains.)
    rp_pidc_lst_f_ = os.path.join('..', 'data', 'NMR', 'multimodel_lists',
                                       'multimod_2713_hetallchains_hom1chain.lst')
    rp_fasta_f_ = os.path.join(_rp_mmseqs_fasta_dir(sub_dir='multimod_2713_hetallchains_hom1chain'),
                               'multimod_2713_hetallchains_hom1chain.fasta')
    pid_stats_pdf, pidc_stats_pdf = generate_stats(sub_dir='multimod_2713_hetallchains_hom1chain',
                               rp_pidc_lst_f=rp_pidc_lst_f_,
                               rp_fasta_f= rp_fasta_f_,
                               run_and_write_mmseqs2=False,
                               run_and_write_rmsd=False, run_and_write_tms=True, use_mmcif=True)

    rp_stats_dst_dir = os.path.join(_rp_nmr_dir(), 'stats', 'multimod_2713_hetallchains_hom1chain')
    os.makedirs(rp_stats_dst_dir, exist_ok=True)
    stats_dst_f = os.path.join(rp_stats_dst_dir, 'mm_2713.csv')
    stats_pdf.to_csv(stats_dst_f, index=False)
