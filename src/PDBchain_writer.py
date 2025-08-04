import os
import MDAnalysis as mda
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import mmseqs2


def write_average_structure(pdb_id: str):
    print(f'pdb_id: {pdb_id}')
    rp_raw_pdbs_dir = mmseqs2.rp_raw_pdbs_dir(het_hom='hethom_combined')
    rp_pdb = os.path.join(rp_raw_pdbs_dir, f'{pdb_id}.pdb')
    rp_pdb_mean_coords_dir = os.path.join(rp_raw_pdbs_dir, 'mean_coords')
    os.makedirs(rp_pdb_mean_coords_dir, exist_ok=True)
    rp_pdb_mean_coords_f = os.path.join(rp_pdb_mean_coords_dir, f'{pdb_id}.pdb')
    pdb_dict = MMCIF2Dict(rp_pdb)



def write_average_structure_using_MDA(pdb_id: str):
    print(f'pdb_id: {pdb_id}')
    rp_raw_pdbs_dir = mmseqs2.rp_raw_pdbs_dir(het_hom='hethom_combined')
    rp_pdb = os.path.join(rp_raw_pdbs_dir, f'{pdb_id}.pdb')
    rp_pdb_mean_coords_dir = os.path.join(rp_raw_pdbs_dir, 'mean_coords')
    os.makedirs(rp_pdb_mean_coords_dir, exist_ok=True)
    rp_pdb_mean_coords_f = os.path.join(rp_pdb_mean_coords_dir, f'{pdb_id}.pdb')

    u = mda.Universe(rp_pdb)

    n_atoms = len(u.atoms)
    n_frames = len(u.trajectory)

    coords_sum = np.zeros((n_atoms, 3))

    for ts in u.trajectory:
        coords_sum += u.atoms.positions

    coords_avg = coords_sum / n_frames
    u.atoms.positions = coords_avg

    u.atoms.write(rp_pdb_mean_coords_f)


if __name__ == '__main__':

    # PATHS OF PDBid_CHAINS:
    rp_mmseqs_results_dir = mmseqs2.rp_mmseqs_results_dir(het_hom='hethom_combined')
    rp_homologues_30_20_90 = os.path.join(rp_mmseqs_results_dir, 'homologues_30_20_90.csv')
    rp_non_homologues_30_20_90 = os.path.join(rp_mmseqs_results_dir, 'non_homologues_30_20_90.csv')

    # PidChains:
    homologues_30_20_90_pdf = pd.read_csv(rp_homologues_30_20_90)
    non_homologues_30_20_90_pdf = pd.read_csv(rp_non_homologues_30_20_90)

    homologues_query = homologues_30_20_90_pdf['query'].tolist()
    homologues_target = homologues_30_20_90_pdf['target'].tolist()
    homologues_all = list(set(homologues_query + homologues_target))
    homologues_all.sort()

    non_homologues_query = non_homologues_30_20_90_pdf['query'].tolist()
    non_homologues_target = non_homologues_30_20_90_pdf['target'].tolist()
    non_homologues_all = list(set(non_homologues_query + non_homologues_target))
    non_homologues_all.sort()

    all_pidChains = list(set(homologues_all + non_homologues_all))
    all_pidChains.sort()

    all_pdbids = list(set([pidchain.split('_')[0] for pidchain in all_pidChains]))
    all_pdbids.sort()

    for pdbid in all_pdbids:
        if pdbid == '1CIR':
            continue
        else:
            write_average_structure_using_MDA(pdbid)

