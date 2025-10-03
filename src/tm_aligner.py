"""
Downloaded via Zhang lab website: http://zhanggroup.org/TM-align/
Compile TMalign_exe.cpp:
 - in src/TMalign_exe/LinuxOS: `$ g++ -O3 -ffast-math -lm -o TMalign_exe TMalign_exe.cpp`
 - in src/TMalign_exe/MacOS which has `#include <malloc.h>` on line 73 commented out.
   Prior to compiling the Mac (Silicon) download, I had to comment out `#include <malloc.h>` on line 10 of `basic_fun.h`
   `// #include <malloc.h> // is replaced in Mac by #include <stdlib.h> which is already on line 6.`
   After commenting it out, I compiled, and finally deleted basic_fun.h and all other files.
"""

import os, glob
from pathlib import Path
from time import time
from tqdm import tqdm
import subprocess
import itertools
import statistics
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
import platform

"""
(CGT4o) TM-score and interpretation:
  1.0 = perfect structural match (identical structures)
> 0.8 = often considered nearly identical structures (small conformational shifts only)
> 0.5 = typically indicates same fold
< 0.2 ≈ random similarity
"""

def _rp_nmr_dir() -> str:
    return os.path.join('..', 'data', 'NMR')

def _rp_tms_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'TMscores', sub_dir)

def _read_all_pdbs_from_txt(txt_f: str) -> list:
    rp_pidchains_lst_f = os.path.join('..', 'data', 'NMR', 'multimodel_lists', txt_f)
    if 'singlemod' == txt_f.split('_')[1]:
        with open(rp_pidchains_lst_f, 'r') as f:
            next(f)  # ignore first line
            pidchains = f.read().split()
    else:
        with open(rp_pidchains_lst_f, 'r') as f:
            pidchains = f.read().split('\n')
            pidchains = pidchains[:-1]
    rp_pdb_files = []
    PDB_dir = os.path.join('..', 'data', 'NMR', 'pdb_chains', 'hethom_combined')  # TODO update
    for pidchain in pidchains:
        rp_pdb_files.append(os.path.join(PDB_dir, f'{pidchain}.pdb'))
    return rp_pdb_files

op_sys = platform.system()  # 'Linux', 'Darwin' (Darwin = macOS), or 'Windows'
TMALIGN_BIN = os.path.join('TMalign_exe', op_sys, 'TMalign_exe')

def compute_tm_from_mp_pool(idx_pair, rp_all_pdb_files):
    """
    Takes a pair of integers, used as indexes to the list of PDB files (given as the relative path strings to the
    appropriate PDB/PDBchain files (which are expected to be single model, or a single mean of all models).)
    Compute TM-score between each pair of PDB/PDBchain files.
    """
    f1_idx, f2_idx = idx_pair
    f1 = rp_all_pdb_files[f1_idx]
    f2 = rp_all_pdb_files[f2_idx]
    f1_pdbname = os.path.basename(f1).removesuffix('.pdb')
    f2_pdbname = os.path.basename(f2).removesuffix('.pdb')

    assert Path(TMALIGN_BIN).exists(), f'{TMALIGN_BIN} does not exist'
    assert os.path.isfile(f1), f'{f1} not found'
    assert os.path.isfile(f2), f'{f2} not found'
    try:
        result = subprocess.run(
            [TMALIGN_BIN, f1, f2],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=60
        )
        output = result.stdout
        for line in output.splitlines():
            if line.startswith('TM-score='):
                tm = float(line.split('=')[1].split()[0])
                print(f'\nTM-score={tm} for {f1_pdbname}:{f2_pdbname}')
                return f1_pdbname, f2_pdbname, tm
    except Exception as e:
        print(f'Error computing TM-score for {f1} vs {f2}: {e}')
        return f1_idx, f2_idx, np.nan

def compute_tm_all_vs_all(rp_pdb_subdir: str):
    pidc = os.path.basename(rp_pdb_subdir).removesuffix('.pdb').split('_')[0]
    print(pidc)
    rp_all_models_pdbs = sorted(glob.glob(os.path.join(rp_pdb_subdir, '*.pdb')))
    N = len(rp_all_models_pdbs)
    # pairs = list(itertools.combinations(range(N), 2))
    # pairs = list(itertools.combinations_with_replacement(range(N), 2))  # including self-self pairs (+ve control).
    pairs = list(itertools.product(range(N), repeat=2))  # including duplicates (0,1) and (1,0).
    print(f'Total pairs to compute: {len(pairs)}')
    n_workers = max(cpu_count() - 1, 1)
    print(f'Using {n_workers} parallel workers.')

    _compute_tm_from_mp_pool = partial(compute_tm_from_mp_pool, rp_all_pdb_files=rp_all_models_pdbs)

    results = {'query': [], 'target': [], 'TMscore': []}
    tms_matrix = np.zeros((N, N))

    # Simple heuristic: bigger groups → bigger chunks, but never below 1
    chunksize = max(1, len(pairs) // (4 * cpu_count()))
    # Optionally clamp to a reasonable upper bound
    chunksize = min(chunksize, 10)

    with Pool(processes=n_workers) as pool:
        for result_ in tqdm(
                pool.imap_unordered(_compute_tm_from_mp_pool, pairs, chunksize=10),
                total=len(pairs), desc="Processing TM-scores"
        ):
            pdb1, pdb2, tms = result_
            if np.isnan(tms):
                continue
            # results['query'].append(pdb1)
            # results['target'].append(pdb2)
            # results['TMscore'].append(tms)
            tms_matrix[i, j] = round(tms, 3)

    # Saving both to dataframe/csv and npz (using the former for readibility, until confident it works as intended):
    rp_dst_tmalign_dir = _rp_tms_dir(sub_dir='multimod_2713_hetallchains_hom1chain')
    rp_dst_tms_mats_dir = os.path.join(rp_dst_tmalign_dir, 'tms_matrices')
    os.makedirs(rp_dst_tms_mats_dir, exist_ok=True)
    # to csv:
    res_pdf = pd.DataFrame(results)
    res_pdf.to_csv(os.path.join(rp_dst_tms_mats_dir, f'{pidc}.csv'), index=False)
    # to npz:
    rp_tms_mat_pidc = os.path.join(rp_dst_tms_mats_dir, pidc)
    np.savez(f'{rp_tms_mat_pidc}.npz', mat=tms_matrix)

def compute_tm(rp_pdb1, rp_pdb2):
    """
    Note: TM-score, unlike my RMSD calculation, does not take in arrays of coordinates.
    Instead, it expects (valid) PDB files.
    Compute TM-score between two pdb files.
    """
    try:
        result = subprocess.run(
            [TMALIGN_BIN, rp_pdb1, rp_pdb2],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=60
        )
        output = result.stdout
        for line in output.splitlines():
            if str(line).startswith('TM-score='):
                tm = float(line.split('=')[1].split()[0])
                return tm
    except Exception as e:
        print(f'Error computing TM-score for {rp_pdb1} vs {rp_pdb2}: {e}')
        return np.nan


if __name__ == '__main__':

    # Get list of pidc to tmalign:
    # dir_ = os.path.join(_rp_tms_dir(sub_dir='multimod_2713_hetallchains_hom1chain'))
    # rp_csv_files = sorted(glob.glob(os.path.join(dir_, '*.csv')))
    # pidc_list = [os.path.basename(rp_csv_f).removesuffix('.csv').removeprefix('TMS_') for rp_csv_f in rp_csv_files]
    # # Make list of pid which include one chain only:
    # from collections import defaultdict
    # pid_c_dict = defaultdict(list)
    # pidc_list.sort()
    # for pidc in pidc_list:
    #     pid, c = pidc.split('_')
    #     pid_c_dict[pid].append(c)
    # pid_c_dict = dict(pid_c_dict)
    #
    # # Get list of pid already done:
    # with open(os.path.join(dir_, 'tms_matrices', 'pid_done.txt'), 'r') as f:
    #     pids_done = f.read().splitlines()
    # pids_done = [pid.removesuffix('.npz') for pid in pids_done]
    #
    # # Get list of pids with only 1 chain (and rename the .npz file to include the chain suffix)
    # only1chain = list()
    # for pid in pids_done:
    #     if pid in pid_c_dict:
    #         if len(pid_c_dict[pid]) == 1:
    #             only1chain.append(pid)
    # print(only1chain)
    #
    # for pid in only1chain:
    #     rp_old_npz = os.path.join(dir_, 'tms_matrices', f'{pid}.npz')
    #     rp_new_npz = os.path.join(dir_, 'tms_matrices', f'{pid}_{pid_c_dict[pid][0]}.npz')
    #     os.rename(rp_old_npz, rp_new_npz)

    # rp_tms_mat = os.path.join(_rp_tms_dir('multimod_2713_hetallchains_hom1chain'), 'tms_matrices', )
    # loaded = np.load(rp_rmsd_mat)
    # rmsd_mat = loaded['mat']
    # dendrgm_heatmap_contourmap(rp_rmsd_mat)


    # (NOTE: I DISCOVERED I COULD NOT RE-USE THE MEAN COORDS (NP ARRAYS) CALCULATED & USED IN COMPUTING RMSDS,
    # BECAUSE TM-ALIGN WANTS PDB FILE INPUTS ONLY.
    # SO, MEAN COORD PDBS FOR TM WERE CALCULATED BY BIOPTYHON VIA src/TMalign_exe/pdb_chain_writer/write_mean_models_to_pdb(),
    # AND SAVED TO data/NMR/pdb_chains/multimod_2713_hetallchains_hom1chain/mean_coords.

    start = time()
    # # CALCULATE TM-SCORES OF EACH MODEL VS THE MEAN OF ALL MODELS:
    rp_pidc_2713_dir = os.path.join(_rp_nmr_dir(), 'pdb_chains', 'multimod_2713_hetallchains_hom1chain', 'per_model')
    rp_subdirs = sorted(glob.glob(os.path.join(rp_pidc_2713_dir, '*')))
    rp_dst_tms_mats_dir = os.path.join(_rp_tms_dir('multimod_2713_hetallchains_hom1chain'), 'tms_matrices')
    os.makedirs(rp_dst_tms_mats_dir, exist_ok=True)

    less_than_3_CA = ['1GAC_A', '1GAC_B', '1GAC_C', '1GAC_D', '1WCO_A',
                      '2AIZ_B', '2K1Q_B', '2M9P_B', '2M9Q_B', '2MX6_B', '3CYS_B']

    rp_subdirs = [rp_subdir for rp_subdir in rp_subdirs if os.path.basename(rp_subdir) not in less_than_3_CA]

    for rp_subdir in rp_subdirs:
        rp_pdbs = sorted(glob.glob(os.path.join(rp_subdir, '*.pdb')))
        tms_matrix = np.zeros((len(rp_pdbs), len(rp_pdbs)))
        pidc = None
        for i, rp_pdb1 in enumerate(rp_pdbs):
            if pidc is None:
                pidc = os.path.basename(rp_pdb1).removesuffix('.pdb')[:-3]
                print(pidc)
            for j, rp_pdb2 in enumerate(rp_pdbs):
                tms_matrix[i][j] = compute_tm(rp_pdb1, rp_pdb2)

        rp_tms_mat_pidc_npz = os.path.join(rp_dst_tms_mats_dir, f'{pidc}.npz')
        np.savez(rp_tms_mat_pidc_npz, mat=tms_matrix)


    # rp_mean_coords_pidchains_dir = os.path.join(rp_pidc_2713_dir, 'mean_coords')
    # rp_mean_coords_pidc_pdbs = sorted(glob.glob(os.path.join(rp_mean_coords_pidchains_dir, '*.pdb')))
    # rp_permodel_pdb_dir = os.path.join(rp_pidc_2713_dir, 'per_model')
    #
    # for rp_mean_coords_pidc_pdb in rp_mean_coords_pidc_pdbs:
    #     tms_per_model = {'pidc': [], 'pidc_model': [], 'TM-score': [], 'min_TMS': [],
    #                      'max_TMS': [], 'mean_TMS': [], 'stdev_TMS': []}
    #     pidc_name = os.path.basename(rp_mean_coords_pidc_pdb).removesuffix('.pdb')
    #     print(pidc_name)
    #     # NOTE: TM-ALIGN RETURNS NONE FOR PROTEINS WITH LESS THAN 3 CA ATOMS:
    #     if (pidc_name == '1GAC_A' or pidc_name == '1GAC_B' or pidc_name == '1GAC_C' or pidc_name == '1GAC_D' or
    #             pidc_name == '1WCO_A' or pidc_name == '2AIZ_B' or pidc_name == '2K1Q_B' or pidc_name == '2M9P_B' or
    #             pidc_name == '2M9Q_B' or pidc_name == '2MX6_B' or pidc_name == '3CYS_B'):
    #         continue
    #     tm_scores = []
    #     rp_permodel_pidc_pdbs = sorted(glob.glob(os.path.join(rp_permodel_pdb_dir, pidc_name, '*.pdb')))
    #     for rp_permodel_pidc_pdb in rp_permodel_pidc_pdbs:
    #         tms = compute_tm(rp_pdb1=rp_permodel_pidc_pdb, rp_pdb2=rp_mean_coords_pidc_pdb)
    #         tms_per_model['pidc_model'].append(os.path.basename(rp_permodel_pidc_pdb).removesuffix('.pdb'))
    #         tm_scores.append(tms)
    #
    #     # Note: I calc all the tms values, assign to dict, which is assigned empty lists for the other 5 values.
    #     # Then I make the pdf. This is to allow me to sort by tms values before adding the scalar values to the top of
    #     # the pdf columns (otherwise if I sorted after adding them, they're end up scatter around different rows
    #     # instead of the top row).
    #     tms_per_model['TM-score'] = tm_scores
    #     empty_list = [np.nan for _ in tm_scores]
    #     tms_per_model['pidc'] = ['' for _ in tm_scores]
    #     tms_per_model['min_TMS'] = empty_list
    #     tms_per_model['max_TMS'] = empty_list
    #     tms_per_model['mean_TMS'] = empty_list
    #     tms_per_model['stdev_TMS'] = empty_list
    #     tms_pdf = pd.DataFrame(tms_per_model)
    #     tms_pdf = tms_pdf.sort_values(by=['TM-score'], ascending=[True])
    #
    #     tms_pdf['pidc'] = [pidc_name] + ['' for _ in tm_scores][:-1]
    #     tms_pdf['min_TMS'] = [round(min(tm_scores), 4)] + empty_list[:-1]
    #     tms_pdf['max_TMS'] = [round(max(tm_scores), 4)] + empty_list[:-1]
    #     tms_pdf['mean_TMS'] = [round(statistics.mean(tm_scores), 4)] + empty_list[:-1]
    #     tms_pdf['stdev_TMS'] = [round(statistics.stdev(tm_scores), 4)] + empty_list[:-1]
    #
    #     rp_pidc_dir = os.path.join(_rp_tmscores_dir(sub_dir='multimod_2713_hetallchains_hom1chain'))
    #     os.makedirs(rp_pidc_dir, exist_ok=True)
    #     tms_pdf.to_csv(os.path.join(rp_pidc_dir, f'TMS_{pidc_name}.csv'), index=False)
    print(f'Completed in {round((time() - start) / 60)} minutes.')
