import os
from time import time
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.Data.IUPACData import protein_letters_3to1
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

AA = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
      "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]


def three_to_one(resname):
    return protein_letters_3to1[resname.capitalize()]


def align_alpha_carbons(model1, model2):
    atom_types = ['CA']
    model1_seq = "".join([three_to_one(r.resname) for r in model1.get_residues() if r.resname in AA])
    model2_seq = "".join([three_to_one(r.resname) for r in model2.get_residues() if r.resname in AA])
    assert len(model1_seq) == len(model2_seq), 'Sequences should be of identical length.'
    model1_coords = [a.coord for a in model1.get_atoms() if a.parent.resname in AA and a.name in atom_types]
    model2_coords = [a.coord for a in model2.get_atoms() if a.parent.resname in AA and a.name in atom_types]
    si = SVDSuperimposer()
    np_model1_coords = np.array(model1_coords)
    np_model2_coords = np.array(model2_coords)
    si.set(np_model1_coords, np_model2_coords)
    si.run()
    return si


def compute_rmsd_matrix(pdbid_chain: str, het_hom: str):
    parser = PDB.MMCIFParser(QUIET=True)
    pdbid, chain = pdbid_chain.split('_')
    structure = parser.get_structure("", f'../data/NMR/raw_cifs/{het_hom}/{pdbid}.cif')

    models = list(structure)
    n_models = len(models)
    print(f'{pdbid_chain} has {n_models} models.')
    rmsd_matrix = np.zeros((n_models, n_models))
    for i in range(n_models):
        for j in range(i + 1, n_models):
            chain_i = models[i][chain]
            chain_j = models[j][chain]
            si = align_alpha_carbons(chain_i, chain_j)
            rmsd = si.get_rms()
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
    print(f'RMSD matrix:\n{rmsd_matrix}')
    return rmsd_matrix, n_models


def build_dendrogram(rmsd_mat, n_models):
    condensed = squareform(rmsd_mat)
    linkage_matrix = linkage(condensed, method='ward')
    plt.figure(figsize=(10, 5))
    dendrogram(linkage_matrix, labels=[f'Model {i}' for i in range(n_models)])
    plt.title('Hierarchical clustering of NMR models (RMSD) for {pdb}')
    plt.ylabel('Distance (Å)')
    plt.tight_layout()
    plt.show()
    return linkage_matrix


def build_heatmap(rmsd_matrix, linkage_matrix):
    sns.clustermap(
        rmsd_matrix,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        cmap="viridis",
        linewidths=0.5,
        figsize=(8, 8),
        cbar_kws={'label': 'RMSD (Å)'}
    )
    plt.title("RMSD heatmap of NMR models")
    plt.show()


def build_contour_map(rmsd_matrix):
    plt.figure(figsize=(8, 6))
    c = plt.contourf(rmsd_matrix, levels=15, cmap="viridis")
    plt.colorbar(c, label="RMSD (Å)")
    plt.xlabel("Model index")
    plt.ylabel("Model index")
    plt.title("RMSD Contour Plot")
    plt.show()


def cluster_models(rmsd_matrix, threshold=2.0):
    """
    Cluster models so any pair differing by more than threshold Å is in different clusters.
    Returns: dict mapping cluster IDs to list of model indices
    """
    condensed = squareform(rmsd_matrix)  # Convert to condensed distance
    Z = linkage(condensed, method='average')  # Hierarchical clustering
    cluster_labels = fcluster(Z, t=threshold, criterion='distance')  # Assign clusters

    clusters = {}
    for i, label in enumerate(cluster_labels):
        clusters.setdefault(label, []).append(i)
    return clusters


def clusters_to_json_dict(clusters):
    """
    Convert clusters dict {id: [indices]} to JSON-friendly dict {ensembleN: [modelN...]}
    """
    json_dict = {}
    for idx, (cluster_id, model_indices) in enumerate(sorted(clusters.items())):
        ensemble_key = f'ensemble{idx + 1}'
        json_dict[ensemble_key] = [f'model{m}' for m in model_indices]
    return json_dict


def write_clusters_json(json_dict, output_path):
    """
    Write the JSON dict to file.
    """
    with open(output_path, 'w') as f:
        json.dump(json_dict, f, indent=4)


def cluster_save_to_json(rmsd_mat):
    clusters = cluster_models(rmsd_mat, threshold=2.0)
    json_dict = clusters_to_json_dict(clusters)
    mm_ensbls = '../data/NMR/multimodel_ensembles'
    os.makedirs(mm_ensbls, exist_ok=True)
    output_file = f'{mm_ensbls}/{pdbid_chain}_ens.json'
    write_clusters_json(json_dict, output_file)
    print(f'Saved ensembles to {output_file}')


def rmsd_reference(np_model1_coords, np_model2_coords) -> float:
    si = SVDSuperimposer()
    si.set(np_model1_coords, np_model2_coords)
    si.run()
    rmsd = si.get_rms()
    return rmsd

# NOTE I DO NOT SUPERIMPOSE COORDINATES OF DIFFERENT MODELS ONTO THE RANDOMLY CHOSEN REFERENCE MODEL,
# PRIOR TO COMPUTING MEAN & STD DEV. IT MIGHT BE BENEFICIAL TO DO THIS BUT IN THIS PRELIMINARY STAGE, I'VE NOT DONE SO.
def mean_structure(pdbid_chain: str):
    rp_pdbid_chain = f'../data/NMR/tokenised_cifs/heteromeric/{pdbid_chain}'
    pdf = pd.read_csv(f'{rp_pdbid_chain}.ssv', sep=' ')
    # pdf = pdf.sort_values(by=["A_pdbx_PDB_model_num", "A_id"])
    model_numbers = pdf['A_pdbx_PDB_model_num'].unique()
    coord_list = []
    ref_A_id = None
    ref_S_mon_id = None
    ref_S_seq_id = None

    for i, model_num in enumerate(model_numbers):
        model_df = pdf[pdf['A_pdbx_PDB_model_num'] == model_num].copy()
        S_mon_id_values = model_df['S_mon_id'].values
        S_seq_id_values = model_df['S_seq_id'].values

        if i == 0: # first model's labels for validation (asserts)
            ref_A_id = model_df['A_id'].values
            ref_S_mon_id = S_mon_id_values
            ref_S_seq_id = S_seq_id_values
        else:
            assert np.array_equal(S_mon_id_values, ref_S_mon_id), f'S_mon_id mismatch in model {model_num}'
            assert np.array_equal(S_seq_id_values, ref_S_seq_id), f'S_seq_id mismatch in model {model_num}'

        coords = model_df[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
        coord_list.append(coords)

    coord_array = np.stack(coord_list)  # Stack into 3D array: (n_models, n_atoms, 3)
    mean_coords = coord_array.mean(axis=0)
    std_coords = coord_array.std(axis=0, ddof=1)

    mean_df = pd.DataFrame({
        'S_mon_id': ref_S_mon_id,
        'S_seq_id': ref_S_seq_id,
        'A_id': ref_A_id,
        'mean_x': mean_coords[:, 0],
        'std_x': std_coords[:, 0],
        'mean_y': mean_coords[:, 1],
        'std_y': std_coords[:, 1],
        'mean_z': mean_coords[:, 2],
        'std_z': std_coords[:, 2]
    })
    cols_to_cast = ['mean_x', 'std_x', 'mean_y', 'std_y', 'mean_z', 'std_z']
    mean_df[cols_to_cast] = mean_df[cols_to_cast].astype(np.float16)  # I assume this loss of precision is ok..
    mean_df = mean_df.reset_index(drop=True)

    rp_rmsd_dir = os.path.join('..', 'data', 'NMR', 'RMSD', 'heteromeric', 'mean_coords')
    os.makedirs(rp_rmsd_dir, exist_ok=True)
    rp_rmsd_csv = os.path.join(rp_rmsd_dir, f'{pdbid_chain}.csv')
    mean_df.to_csv(rp_rmsd_csv, index=False)
    return mean_df


if __name__ == '__main__':
    _meric = 'heteromeric'

    with open(f'../data/NMR/multimodel_lists/{_meric[:3]}_multimod_2104_pdbid_chains.txt', 'r') as f:
        pdbid_chains = f.readlines()

    pdbid_chain = '1A0N_A'

    rmsd_mat, n_models = compute_rmsd_matrix(pdbid_chain, het_hom=_meric)
    # # linkage_matrix = build_dendrogram(rmsd_mat, n_models)
    # # build_heatmap(rmsd_mat, linkage_matrix)
    # # build_contour_map(rmsd_mat)

    # cluster_save_to_json(rmsd_mat)
    #
    # # clusters = cluster_models(rmsd_mat, threshold=2.0)
    # # json_dict = clusters_to_json_dict(clusters)
    # # mm_ensbls = '../data/NMR/multimodel_ensembles'
    # # os.makedirs(mm_ensbls, exist_ok=True)
    # # output_file = f'{mm_ensbls}/{pdbid_chain}_ens.json'
    # # write_clusters_json(json_dict, output_file)
    # # print(f'Saved ensembles to {output_file}')

    df = mean_structure(pdbid_chain)