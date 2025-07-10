import os
import json
import numpy as np
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
    assert len(model1_seq) == len(model2_seq), "Sequences should be of identical length."
    model1_coords = [a.coord for a in model1.get_atoms() if a.parent.resname in AA and a.name in atom_types]
    model2_coords = [a.coord for a in model2.get_atoms() if a.parent.resname in AA and a.name in atom_types]
    si = SVDSuperimposer()
    si.set(np.array(model1_coords), np.array(model2_coords))
    si.run()
    return si


def compute_rmsd_matrix(pdbid_chain):
    parser = PDB.MMCIFParser(QUIET=True)
    pdbid, chain = pdbid_chain.split('_')
    structure = parser.get_structure("", f'../data/NMR/raw_cifs/heteromeric/{pdbid}.cif')

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
    # Convert to condensed distance:
    condensed = squareform(rmsd_matrix)
    # Hierarchical clustering:
    Z = linkage(condensed, method='average')
    # Assign clusters:
    cluster_labels = fcluster(Z, t=threshold, criterion='distance')
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


if __name__ == "__main__":
    with open('../data/NMR/multimodel_PDBid_chains/het_multimod_2104_pdbid_chains.txt', 'r') as f:
        pdbid_chains = f.readlines()

    pdbid_chain = '1A0N_A'
    rmsd_mat, n_models = compute_rmsd_matrix(pdbid_chain)
    # linkage_matrix = build_dendrogram(rmsd_mat, n_models)
    # build_heatmap(rmsd_mat, linkage_matrix)
    # build_contour_map(rmsd_mat)
    clusters = cluster_models(rmsd_mat, threshold=2.0)
    json_dict = clusters_to_json_dict(clusters)
    dst_dir = '../data/NMR/multimodel_PDBid_chains/ensembles'
    os.makedirs(dst_dir, exist_ok=True)
    output_file = f'{dst_dir}/{pdbid_chain}_ens.json'
    write_clusters_json(json_dict, output_file)
    print(f'Saved ensembles to {output_file}')
