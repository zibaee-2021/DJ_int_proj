import os, glob
import statistics as stats
from time import time
from typing import Tuple
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
import RMSD
import mmseqs2
import tm_aligner

# RELATIVE PATHS:
def _rp_nmr_dir() -> str:
    return os.path.join('..', 'data', 'NMR')

def _rp_raw_cifs_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'raw_cifs', sub_dir)

def _rp_stats_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'stats', sub_dir)

def _rp_rmsd_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'RMSD', sub_dir)

def _rp_pdb_chains_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'pdb_chains', sub_dir)

def _rp_tmscores_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'TMscores', sub_dir)

def _rp_parsed_cifs_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', sub_dir)

def _rp_mmseqs_dir(sub_dir) -> str:
    return os.path.join(_rp_nmr_dir(), 'mmseqs', sub_dir)

def _rp_mmseqs_fasta_dir(sub_dir) -> str:
    return os.path.join(_rp_mmseqs_dir(sub_dir), 'fasta')

def _calc_identity_for_stats(het_hom: str, pdbid: str, filt_pdf) -> str:
    if not filt_pdf.empty:
        print(f'There is some heterogeneity of sequence between some chains. '
              f'See data/NMR/mmseqs/{het_hom}/results/{pdbid}.csv')
        identity = '(chains not identical, see mmseqs pdf)'
    else:
        print('All chains have high identity, if not 100%')
        identity = '0'
    return identity

# TEMPORARY DUPLICATE OF THAT IN RMSD.py
def _write_rmsds_and_stats(sub_dir: str):
    rp_mean_coords_pidc_dir = os.path.join(_rp_rmsd_dir(sub_dir), 'mean_coords')
    rp_mean_coords_pidc_csvs = sorted(glob.glob(os.path.join(rp_mean_coords_pidc_dir, '*.csv')))

    for rp_mean_coords_pidc_csv in rp_mean_coords_pidc_csvs:
        rmsd_per_model = {'pidc': [], 'pidc_model': [], 'RMSD': [],
                          'min_RMSD': [], 'max_RMSD': [], 'mean_RMSD': [], 'stdev_RMSD': []}
        pidc = os.path.basename(rp_mean_coords_pidc_csv).removesuffix('.csv')
        rmsd_per_model['pidc'].append(pidc)
        rmsds = []
        mean_pidc_pdf = pd.read_csv(rp_mean_coords_pidc_csv)
        mean_coords_pidc = mean_pidc_pdf[['mean_x', 'mean_y', 'mean_z']].values
        rp_parsed_cif_ssv = os.path.join(_rp_parsed_cifs_dir(sub_dir), f'{pidc}.ssv')
        pidc_pdf = pd.read_csv(os.path.join(rp_parsed_cif_ssv), sep=' ')
        for model_num, pidc_model in pidc_pdf.groupby('A_pdbx_PDB_model_num'):
            pidc_coords = pidc_model[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
            if pidc == '1MSH_A' and model_num == 30:
                continue
            elif pidc == '2JSC_A' and model_num == 1:
                continue
            else:
                rmsd = RMSD.calculate_rmsd(coords1=mean_coords_pidc, coords2=pidc_coords)
                rmsds.append(round(rmsd, 4))
                pidc_model = f'{pidc}_{model_num}' if model_num >= 10 else f'{pidc}_0{model_num}'
                rmsd_per_model['pidc_model'].append(pidc_model)

        if len(rmsds) == 0:
            print(f'rmsds for {pidc} is empty: {rmsds}. Assigning np.nan.')
            min_rmsd, max_rmsd = np.nan, np.nan
            mean_rmsd, stdev_rmsd = np.nan, np.nan
        else:
            min_rmsd, max_rmsd = min(rmsds), max(rmsds)
            mean_rmsd, stdev_rmsd = stats.mean(rmsds), stats.stdev(rmsds)

        # Note: I build the rmsds, but otherwise empty, dict & pdf first, to be able to sort by this before adding
        # the remaining, whcih are scalar or single string values, which are all on just the top line, otherwise
        # it'd look weird with these single row values randomly on a different line, with nans everywhere else:
        rmsd_per_model['RMSD'] = rmsds
        empty_list = [np.nan for _ in rmsds]
        rmsd_per_model['pidc'] = ['' for _ in rmsds]
        rmsd_per_model['min_RMSD'] = empty_list
        rmsd_per_model['max_RMSD'] = empty_list
        rmsd_per_model['mean_RMSD'] = empty_list
        rmsd_per_model['stdev_RMSD'] = empty_list
        rmsd_pdf = pd.DataFrame(rmsd_per_model)
        rmsd_pdf = rmsd_pdf.sort_values(by=['RMSD'], ascending=[False])

        rmsd_pdf['pidc'] = [pidc] + ['' for _ in rmsds][: -1]
        rmsd_pdf['min_RMSD'] = [round(min_rmsd, 4)] + empty_list[: -1]
        rmsd_pdf['max_RMSD'] = [round(max_rmsd, 4)] + empty_list[: -1]
        rmsd_pdf['mean_RMSD'] = [round(mean_rmsd, 4)] + empty_list[: -1]
        rmsd_pdf['stdev_RMSD'] = [round(stdev_rmsd, 4)] + empty_list[: -1]

        rp_pidc_dir = os.path.join(_rp_rmsd_dir('multimod_2713_hetallchains_hom1chain'))
        os.makedirs(rp_pidc_dir, exist_ok=True)
        rmsd_pdf.to_csv(os.path.join(rp_pidc_dir, f'RMSD_{pidc}.csv'), index=False)

def _compute_tms(run_and_write_tms: bool, pidc: str, sub_dir: str) -> tuple:

    if not run_and_write_tms: # THEN READ IN ALREADY CALCULATED VALUES
        # NOTE: TM-ALIGN RETURNS NONE FOR PROTEINS WITH LESS THAN 3 CA ATOMS (WHICH THESE 11 PDB-CHAINS HAVE):
        if (pidc == '1GAC_A' or pidc == '1GAC_B' or pidc == '1GAC_C' or pidc == '1GAC_D' or
                pidc == '1WCO_A' or pidc == '2AIZ_B' or pidc == '2K1Q_B' or pidc == '2M9P_B' or
                pidc == '2M9Q_B' or pidc == '2MX6_B' or pidc == '3CYS_B'):
            min_tms, max_tms, mean_tms, stdev_tms = np.nan, np.nan, np.nan, np.nan
        else:
            rp_tmscores_dir = tm_aligner._rp_tmscores_dir(sub_dir)
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
        else:
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

def _compute_rmsd(run_and_write_rmsd: bool, pidc: str, sub_dir: str):
    if run_and_write_rmsd:
        _write_rmsds_and_stats(sub_dir)
    rp_rmsd_per_model_dir = os.path.join(_rp_rmsd_dir(sub_dir))
    rp_rmsd_per_model_csv_f = os.path.join(rp_rmsd_per_model_dir, f'RMSD_{pidc}.csv')
    rmsd_pdf = pd.read_csv(rp_rmsd_per_model_csv_f)
    rmsds = rmsd_pdf['RMSD'].values
    if rmsds.size > 1:
        min_rmsd, max_rmsd = round(np.min(rmsds), 4), round(np.max(rmsds), 4)
        mean_rmsd, stdev_rmsd = round(stats.mean(rmsds), 4), round(stats.stdev(rmsds), 4)
    else:
        min_rmsd, max_rmsd, mean_rmsd, stdev_rmsd = np.nan, np.nan, np.nan, np.nan
    return min_rmsd, max_rmsd, mean_rmsd, stdev_rmsd

def _allchains_year_and_head(rp_raw_struct_f: str, biopython_parser) -> Tuple[list, int, str]:
    """
    # The following code here in the docstring is not used, for the reason given below:
    # cif_dict = MMCIF2Dict.MMCIF2Dict(rpath_cif)
    # year = cif_dict.get('_citation.year', ['NA'])[0]
    # A few of these return "?", so I chose to use 'MMCIFParser().get_structure(cif)' instead,
    # which often gives same year but sometimes differs by as much as ~3 years.
    """
    # Note: Parser.get_structure() is quite slow, 15-20 secs:
    start = time()
    bio_struct = biopython_parser.get_structure('', rp_raw_struct_f)
    allchains = list(next(bio_struct.get_models()).get_chains())
    all_chains = []
    for ch in allchains:
        all_chains.append(ch.id)
    year = bio_struct.header['deposition_date'][:4]
    head = bio_struct.header['head']
    secs_taken = round(time() - start)
    if secs_taken > 10:
        print(f' bio_parser.get_structure() took {secs_taken} seconds')
    return all_chains, int(year), head

def _pid2c(pidc_list: list) -> dict:
    """
    Takes relative path of txt file that has list of sol NMR PDBid_chains with > 1 model.
    Makes a dict with unique PDBid as the key, mapped to all its chains, according to the list in the txt file.
    Returns a list of the PDBid-chains, and a dictionary mapping PDB ids to chains.
    """
    pid_c_dict = defaultdict(list)
    pidc_list.sort()
    for pidc in pidc_list:
        pid, c = pidc.split('_')
        pid_c_dict[pid].append(c)
    return dict(pid_c_dict)

def _read_multimodel_pidc(rp_pidc: str) -> list:
    if rp_pidc.endswith('.txt'):
        with open(rp_pidc, 'r') as f:
            pidc_list = f.read().split()
    else:  # .lst file
        with open(rp_pidc, 'r') as f:
            pidc_list = f.read().splitlines()
    pidc_list.sort()
    return pidc_list

def _all_prot_only_chains() -> dict:
    """
    Relies on the fact that the lists in the text files hard-coded in this function, are only those that have come
    through the parsing of raw cifs first, which will have removed any PDB chains that have no natural residue
    (alpha-carbons).
    Returns a dict of PDB ids mapped to the number of chains.
    """
    # ONE PIDC PER LINE:
    def _read_mm_pidc(rp_f: str, expected: int) -> list:
        with open(rp_f, 'r') as f:
            lines = f.read().split()
        lines = [line.removesuffix('\n') for line in lines]
        assert len(lines) == expected, f'{len(lines)} is not the expected number ({expected}) of pidc.'
        return lines

    rp_txt_f = os.path.join(_rp_nmr_dir(), 'multimodel_lists', 'het_multimod_2104_PidChains.txt')
    het_pidc = _read_mm_pidc(rp_txt_f, expected=2104)
    het_pid = list(set([het_p[:-2] for het_p in het_pidc]))
    rp_txt_f = os.path.join(_rp_nmr_dir(), 'multimodel_lists', 'hom_multimod_1421_PidChains.txt')
    hom_pidc = _read_mm_pidc(rp_txt_f, expected=1421)
    hom_pid = list(set([hom_p[:-2] for hom_p in hom_pidc]))
    # all_pid = sorted(list(set(het_pid + hom_pid))) # 1541 multimodel pid
    all_pidc = sorted(list(set(het_pidc + hom_pidc))) # 3525 multimodel pidc
    return _pid2c(all_pidc) # 1541 keys (pid)
    # pid_c_count_dict = defaultdict(int)
    # for pid, c in pidc_dict.items():
    #     pid_c_count_dict[pid] = len(pidc_dict[pid])
    # return dict(pid_c_count_dict)

def generate_stats(sub_dir: str, rp_pidc_lst_f: str, rp_fasta_f: str, _chains_year_head_prewritten=True,
                   run_and_write_mmseqs2=False, run_and_write_rmsd=False, run_and_write_tms=False, use_mmcif=True) -> tuple:
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
    pid_stats = {'Pid': [], 'head': [], 'all_chains': [], 'hethom': [], 'year': [], 'all_prot_only_chains': []}
    pid_list, head_list, all_chains_list, hethom_list, year_list, all_prot_only_chains_list= [], [], [], [], [], []

    pidc_stats = {'Pidc': [], 'CA_count': [], 'model_count': [],
                  'minRMSD': [], 'maxRMSD': [], 'meanRMSD': [], 'stdevRMSD': [],
                  'minTMS': [], 'maxTMS': [], 'meanTMS': [], 'stdevTMS': []}#, 'homologues': []}
    pidc_list, ca_count_list, model_count_list = [], [], []
    min_rmsd_list, max_rmsd_list, mean_rmsd_list, stdev_rmsd_list = [], [], [], []
    min_tms_list, max_tms_list, mean_tms_list, stdev_tms_list = [], [], [], []
    year_list, homologues_list = [], []

    if use_mmcif:
        pdb_or_cif = 'cif'
        parser = MMCIFParser(QUIET=True, auth_chains=False)  # Used in _total_chain_count_and_year() below (~line 150).
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
    rp_mmseqs2_results_dir = os.path.join(_rp_mmseqs_dir(sub_dir), 'results')
    os.makedirs(rp_mmseqs2_results_dir, exist_ok=True)
    fname = os.path.basename(rp_fasta_f).removesuffix('.fasta')
    rp_mmseqs2_csv_f = os.path.join(rp_mmseqs2_results_dir, f'{fname}.csv')

    if run_and_write_mmseqs2:
        mmseq2_pdf = mmseqs2.run_mmseqs2_esysrch_cmd_all_vs_all(rp_fasta_f)
        mmseq2_pdf = mmseqs2.dedupe_rm_self_aligns(mmseq2_pdf)
        mmseq2_pdf.to_csv(rp_mmseqs2_csv_f, index=False)  # shape=(33074, 7)
    else:
        mmseq2_pdf = pd.read_csv(rp_mmseqs2_csv_f)  # shape=(33074, 7)

    all_prot_only_chains_dict = _all_prot_only_chains()

    pid_chains_dict = _pid2c(_read_multimodel_pidc(rp_pidc_lst_f))
    for i, (pid, chains) in enumerate(pid_chains_dict.items()):
        if use_mmcif:
            assert pid in raw_pids, f'{pid}.cif not round in raw_cifs.'
        else:
            if pid not in ['9D9B', '7ZE0', '9D9C', '9D9A']:  # (These 4 are not found in 'legacy' PDB files on RCSB)
                assert pid in raw_pids, f'{pid}.pdb not round in raw_pdbs.'
            else:
                print(f'{pid}.pdb is not available in legacy PDB for downloading from RCSB. So, cannot include it.')
                continue

        het_hom = 'het' if pid in raw_het_pids else 'hom'
        rp_raw_struct_dir = rp_raw_struct_hom_dir if het_hom == 'hom' else rp_raw_struct_het_dir
        rp_raw_struct_f = os.path.join(rp_raw_struct_dir, f'{pid}.{pdb_or_cif}')

        allchains, year, head = _allchains_year_and_head(rp_raw_struct_f, parser)
        all_chains_list.append(allchains)
        year_list.append(year)
        head_list.append(head)

        all_prot_only_chains = all_prot_only_chains_dict[pid]
        all_prot_only_chains_list.append(all_prot_only_chains)
        pid_list.append(pid)
        hethom_list.append(het_hom)

        for chain in chains:
            pidc = f'{pid}_{chain}'
            print(pidc)
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
            rp_parsed_cif_ssv_f = os.path.join(_rp_parsed_cifs_dir(sub_dir), f'{pidc}.ssv')
            pidc_pdf = pd.read_csv(rp_parsed_cif_ssv_f, sep=' ')
            model_count = len(pidc_pdf['A_pdbx_PDB_model_num'].unique())
            model_count_list.append(model_count)
            ca_count_list.append(int(pidc_pdf.shape[0] / model_count))

            # GET LIST OF HOMOLOGOUS PDB-CHAHINS FOR GIVEN PDB-CHAIN:
            homologues = mmseqs2.find_homologues_30_20_90(mmseq2_pdf, pidc)
            homologues_list.append(homologues)

        assert (len(pidc_list) == len(ca_count_list) == len(model_count_list) == len(min_rmsd_list) == len(max_rmsd_list) == len(mean_rmsd_list)
                == len(stdev_rmsd_list) == len(min_tms_list) == len(max_tms_list) == len(mean_tms_list)
                == len(stdev_tms_list) == len(homologues_list)), \
            'lengths of column-bound lists for pidc_pdf are not equal.'
        pidc_stats['Pidc'] = pidc_list
        pidc_stats['CA_count'] = ca_count_list
        pidc_stats['model_count'] = model_count_list
        pidc_stats['minRMSD'] = min_rmsd_list
        pidc_stats['maxRMSD'] = max_rmsd_list
        pidc_stats['meanRMSD'] = mean_rmsd_list
        pidc_stats['stdevRMSD'] = stdev_rmsd_list
        pidc_stats['minTMS'] = min_tms_list
        pidc_stats['maxTMS'] = max_tms_list
        pidc_stats['meanTMS'] = mean_tms_list
        pidc_stats['stdevTMS'] = stdev_tms_list
        pidc_stats['homologues'] = homologues_list

    pidc_pdf = pd.DataFrame(pidc_stats)
    pidc_pdf = pidc_pdf.sort_values(by=['maxRMSD', 'maxTMS'], ascending=[False, True])

    pid_stats['Pid'] = pid_list
    pid_stats['head'] = head_list
    pid_stats['all_chains'] = all_chains_list
    pid_stats['hethom'] = hethom_list
    pid_stats['year'] = year_list
    pid_stats['all_prot_only_chains'] = all_prot_only_chains_list
    assert len(pid_list) == len(all_chains_list) == len(hethom_list) == len(year_list), \
        'lengths of column-bound lists for pid_pdf are not equal.'
    pid_pdf = pd.DataFrame(pid_stats)
    pid_pdf = pid_pdf.sort_values(by=['year'], ascending=[True])
    print(f'Completed {len(pidc_list)} PDB-chains in {round((time() - start) / 60)} mins')  # 2713 pidc in 50 mins.
    return pid_pdf, pidc_pdf



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



def _tabulate_year_chain_model_counts(pidc_list: list):
    start = time()
    pid_list, head_list, all_chains_list, hethom_list, year_list, model_count_list = [], [], [], [], [], []
    result_dict = {'pid':pid_list, 'head': head_list, 'all_chains': all_chains_list, 'hethom': hethom_list,
              'year': year_list, 'model_count': model_count_list}

    rp_raw_het_cifs = glob.glob(os.path.join(_rp_raw_cifs_dir('heteromeric'), f'*.cif'))
    het_pids = sorted([os.path.basename(rp_het_cif).removesuffix(f'.cif') for rp_het_cif in rp_raw_het_cifs])
    biopy_parser = MMCIFParser(QUIET=True, auth_chains=False)

    _pid_list = sorted(list(set([pidc[:-2] for pidc in pidc_list])))
    for pid in _pid_list:
        print(pid)
        pid_list.append(pid)
        het_hom = 'het' if pid in het_pids else 'hom'
        hethom_list.append(het_hom)
        subdir = 'heteromeric' if het_hom == 'het' else 'homomeric'
        rp_raw_cif_f = os.path.join(_rp_raw_cifs_dir(subdir), f'{pid}.cif')
        bio_struct = biopy_parser.get_structure('', rp_raw_cif_f)
        model_count = len(list(bio_struct.get_models()))
        model_count_list.append(model_count)
        all_chains_ = list(next(bio_struct.get_models()).get_chains())
        all_chains = []
        for ch in all_chains_:
            all_chains.append(ch.id)
        all_chains_list.append(all_chains)
        head = bio_struct.header['head']
        head_list.append(head)
        year = bio_struct.header['deposition_date'][:4]
        year_list.append(year)

    result_dict['pid'] = pid_list
    result_dict['head'] = head_list
    result_dict['all_chains'] = all_chains_list
    result_dict['hethom'] = hethom_list
    result_dict['year'] = year_list
    result_dict['model_count'] = model_count_list

    result_pdf = pd.DataFrame(result_dict)
    result_pdf.to_csv(os.path.join(_rp_stats_dir('multimod_2713_hetallchains_hom1chain'),
                                   'allchains_head_year_model_hethom.csv'), index=False)
    print(f'Extracted data and written to csv for {len(pid_list)} PDBs in {round((time() - start) / 60)} mins')



if __name__ == '__main__':
    # pid_c_count_dict = _total_natural_prot_chain_count()


    # rp_raw_cifs = sorted(glob.glob(os.path.join(_rp_nmr_dir(), 'raw_cifs', 'hethom_combined', '*.cif')))

    # # VISUALISATIONS:

    # # 1. FASTA SEQUENCE DISTRIBUTION:
    # rp_fasta_f_ = os.path.join(_rp_mmseqs_fasta_dir(sub_dir='multimod_2713_hetallchains_hom1chain'),
    #                            'multimod_2713_hetallchains_hom1chain.fasta')
    # plot_fasta_size_distribution(rp_fasta_f_, x_limit_220=False)

    # # 2. BAR CHART OF ALPHA-CARBON COUNT FOR EACH PDB:
    # # 2.A. CALCULATE CA COUNT FROM PARSED CIFS DATA AND WRITE TO STATS/...LST FILE:
    # rp_parsed_cifs_ssvs_ = sorted(glob.glob(os.path.join(_rp_parsed_cifs_dir('multimod_2713_hetallchains_hom1chain'),
    #                                                      '*.ssv')))
    # ca_counts_ = _calc_ca_counts(rp_parsed_cifs_ssvs_)
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
    # rp_parsed_cifs_ssvs_ = sorted(glob.glob(os.path.join(_rp_parsed_cifs_dir('multimod_2713_hetallchains_hom1chain'),
    #                                                      '*.ssv')))
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

    # # 4. WRITE YEAR (DEPOSITION), MODEL COUNT, CHAIN COUNT - FROM ALL DOWNLOADED NMR CIF FILES (NOT THE PDBCHAINS):
    # # 4.A. EXTRACT YEAR, MODEL COUNT, CHAIN COUNT, FROM PDB FILE AND WRITE TO STATS/...CSV:
    # rp_pidchains_lst_f_ = os.path.join(_rp_nmr_dir(), 'multimodel_lists',
    #                                    'multimod_2713_hetallchains_hom1chain.lst')
    # with open(rp_pidchains_lst_f_, 'r') as f:
    #     pidchains_2713 = sorted(f.read().splitlines())
    # _tabulate_year_chain_model_counts(pidchains_2713)


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


    # I have no clue why '7CLV_A', '7CLV_B', '8J4I_A', '6UJV_A' gave 'null' rmsd stats values from generate_stats() function
    # output to mm_2713_pidc.csv, because they all return values from the _compute_rmsd() function below:
    # min_rmsd, max_rmsd, mean_rmsd, stdev_rmsd = _compute_rmsd(run_and_write_rmsd=True, pidc='6UJV_A',
    #                                                           sub_dir='multimod_2713_hetallchains_hom1chain')
    # pass

    # # MAIN STATS FUNCTION  #########################################################################################
    # # GENERATE STATS PDF AND WRITE TO CSV: (Takes 18 mins to complete 2713 PDB-chains.)
    mm_2713 = 'multimod_2713_hetallchains_hom1chain'
    rp_pidc_lst_f_ = os.path.join('..', 'data', 'NMR', 'multimodel_lists', f'{mm_2713}.lst')
    rp_fasta_f_ = os.path.join(_rp_mmseqs_fasta_dir(sub_dir=mm_2713), f'{mm_2713}.fasta')
    pid_stats_pdf, pidc_stats_pdf = generate_stats(sub_dir=mm_2713, rp_pidc_lst_f=rp_pidc_lst_f_,
                                                   rp_fasta_f= rp_fasta_f_, _chains_year_head_prewritten=False,
                                                   run_and_write_mmseqs2=False, run_and_write_rmsd=False,
                                                   run_and_write_tms=False, use_mmcif=True)

    rp_stats_dst_dir = os.path.join(_rp_nmr_dir(), 'stats', 'multimod_2713_hetallchains_hom1chain')
    os.makedirs(rp_stats_dst_dir, exist_ok=True)

    rp_pid_stats_dst_f = os.path.join(rp_stats_dst_dir, 'mm_2713_pid.csv')
    pid_stats_pdf.to_csv(rp_pid_stats_dst_f, index=False)

    rp_pidc_stats_dst_f = os.path.join(rp_stats_dst_dir, 'mm_2713_pidc.csv')
    pidc_stats_pdf.to_csv(rp_pidc_stats_dst_f, index=False)
    #################################################################################################################