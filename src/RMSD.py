import os, glob
from asyncio import PidfdChildWatcher
from time import time
import statistics
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import mmseqs2
from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.Data.IUPACData import protein_letters_3to1
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

"""(CGT5)

RMSD (Å) and interpretation:

< 1.0 Å	    = Essentially identical; only minor coordinate fluctuations
1.0 – 2.0 Å	= Very similar; typically same fold, possibly different side-chain conformations
2.0 – 4.0 Å	= Similar global fold; might differ in loops or flexible regions
4.0 – 6.0 Å	= Moderate similarity; possibly different topologies or domain orientations
> 6.0 Å	    = Likely different structures; may be different folds entirely
> 10 Å	    = Almost certainly structurally unrelated
"""


AA = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
      "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]


def three_to_one(resname):
    return protein_letters_3to1[resname.capitalize()]

# BUILDING RELATIVE PATHS:
def _rp_nmr_dir():
    return os.path.join('..', 'data', 'NMR')

def _rp_rmsd_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'RMSD', sub_dir)

def _rp_parsed_cifs_dir(sub_dir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', sub_dir)

def rp_mmseqs_dir(het_hom) -> str:
    return os.path.join(_rp_nmr_dir(), 'mmseqs', het_hom)

def rp_mmseqs_results_dir(het_hom) -> str:
    return os.path.join(rp_mmseqs_dir(het_hom), 'results')


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
    structure = parser.get_structure('', f'../data/NMR/raw_cifs/{het_hom}/{pdbid}.cif')

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


def cluster_save_to_json(rmsd_mat, pdbid_chain):
    clusters = cluster_models(rmsd_mat, threshold=2.0)
    json_dict = clusters_to_json_dict(clusters)
    mm_ensbls = '../data/NMR/multimodel_ensembles'
    os.makedirs(mm_ensbls, exist_ok=True)
    output_file = f'{mm_ensbls}/{pdbid_chain}_ens.json'
    write_clusters_json(json_dict, output_file)
    print(f'Saved ensembles to {output_file}')


def calculate_rmsd(coords1, coords2) -> float:
    """
    Calculate RMSD between the two given coordinates (as numpy arrays).
    Uses SVDSuperimposer(). This calculates optimal least-squares superposition (Kabsch algorithm), which is the
    Euclidean RMSD between matched atom pairs. The units of the calculated RMSD value are the same as the units of the
    input coordinates, i.e. Angstroms.
    Lower RMSDs indicate high structural similarity.

    WARNING: The following was from CGT4o:
    < 1 Angs = essentially identical, only minor coordinate fluctuations.
    1-2 Angs = very similar, typically same fold, possibly different side-chain conformations.
    2-4 Angs = Similar global fold; might differ in loops or flexible regions.
    4-6 Angs = Moderate similarity; possibly different topologies or domain orientations.
    > 6 Angs = Likely different structures; may be different folds entirely.
    > 10 Angs = Almost certainly structurally unrelated.
    """
    si = SVDSuperimposer()
    si.set(coords1, coords2)
    si.run()
    rmsd = si.get_rms()
    return rmsd


def calc_rmsds_pidchain_pairs(pidchain1_name: str, pidchain2_name: str, rp_rmsd_mean_coords_dir: str):
    # rmsd_mat, n_models = RMSD.compute_rmsd_matrix(pdbid_chain, het_hom)
    # clusters = RMSD.cluster_models(rmsd_mat, threshold=2.0)
    print(f'pidChain1={pidchain1_name}, pidChain2={pidchain2_name}')
    rp_pidchain1_ssv = os.path.join(rp_rmsd_mean_coords_dir, f'{pidchain1_name}.csv')
    pdbid_chain1_pdf = pd.read_csv(rp_pidchain1_ssv)
    rp_pidchain2_ssv = os.path.join(rp_rmsd_mean_coords_dir, f'{pidchain2_name}.csv')
    pdbid_chain2_pdf = pd.read_csv(rp_pidchain2_ssv)

    coords1 = pdbid_chain1_pdf[['mean_x', 'mean_y', 'mean_z']].values
    coords2 = pdbid_chain2_pdf[['mean_x', 'mean_y', 'mean_z']].values

    rmsd_result = calculate_rmsd(coords1, coords2)
    return rmsd_result



def calc_rmsds_of_models(rp_parsed_cifs_ssv: str, rp_mean_coords_csv: str) -> tuple:
    """
    Calculate RMSD of each model the given PDBchain, vs the mean coordinates of this PDBchain.
    """
    # rmsd_mat, n_models = RMSD.compute_rmsd_matrix(pidChain, het_hom)
    # clusters = RMSD.cluster_models(rmsd_mat, threshold=2.0)
    # ref_structure = mean_stddev_struct_(het_hom, pidChain)
    pidchain_name = os.path.basename(rp_mean_coords_csv).removesuffix('.csv')
    print(pidchain_name)
    ref_structure_pdf = pd.read_csv(rp_mean_coords_csv)
    reference_coords = ref_structure_pdf[['mean_x', 'mean_y', 'mean_z']].values
    pdbid_chain_pdf = pd.read_csv(rp_parsed_cifs_ssv, sep=' ')
    rmsds, model_nums = [], []
    pdc4pdf = [pidchain_name]  # Just for aesthetics of tabular RMSD results.
    for model_num, pdbid_chain_model in pdbid_chain_pdf.groupby('A_pdbx_PDB_model_num'):
        pidchain_coords = pdbid_chain_model[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
        if pidchain_name == '1MSH_A' and model_num == 30:
            continue
        elif pidchain_name == '2JSC_A' and model_num == 1:
            continue
        elif (pidchain_name == '6UJV_A' or pidchain_name == '6UJV_B' or
              pidchain_name == '6UJV_C' or pidchain_name == '7CLV_A' or
              pidchain_name == '7CLV_B' or pidchain_name == '8J4I_A' or
              pidchain_name == '8J4I_B' or pidchain_name == '8J4I_C' or
              pidchain_name == '8J4I_D' or pidchain_name == '8J4I_E' or
              pidchain_name == '8J4I_F'):
            continue
        else:
            rmsd_pidchain = calculate_rmsd(coords1=reference_coords, coords2=pidchain_coords)
            model_nums.append(model_num)
            rmsds.append(rmsd_pidchain)
            pdc4pdf.append('')  # so that the resulting stats table has the pdbid just on the first row only.
    pdc4pdf = pdc4pdf[:-1]
    rmsds = np.array(rmsds, dtype=np.float16)
    return rmsds, model_nums, pdc4pdf

# COPIED OVER FROM pdb_model_stats.py and ammended to match tm_aligner output format:
def _calc_rmsds_and_stats():
    """
    TODO Note that rmsds for 6UJV_A, 7CLV_A, 7CLV_B & 8J4I_A are empty.. need to have a closer look at this to see why...
    """
    start = time()
    rp_mean_coords_pidchains_dir = os.path.join(_rp_rmsd_dir('multimod_2713_hetallchains_hom1chain'), 'mean_coords')
    rp_mean_coords_pidc_csvs = sorted(glob.glob(os.path.join(rp_mean_coords_pidchains_dir, '*.csv')))

    for rp_mean_coords_pidc_csv in rp_mean_coords_pidc_csvs:
        rmsd_per_model = {'pidc': [], 'pidc_model': [], 'RMSD': [],
                          'min_RMSD': [], 'max_RMSD': [], 'mean_RMSD': [], 'stdev_RMSD': []}
        pidc_name = os.path.basename(rp_mean_coords_pidc_csv).removesuffix('.csv')
        print(pidc_name)
        rmsd_per_model['pidc'].append(pidc_name)
        rmsds = []
        mean_pidc_pdf = pd.read_csv(rp_mean_coords_pidc_csv)
        mean_coords_pidc = mean_pidc_pdf[['mean_x', 'mean_y', 'mean_z']].values
        rp_parsed_cif_ssv = os.path.join(_rp_parsed_cifs_dir('multimod_2713_hetallchains_hom1chain'),
                                         f'{pidc_name}.ssv')
        pidc_pdf = pd.read_csv(os.path.join(rp_parsed_cif_ssv), sep=' ')
        for model_num, pdbid_chain_model in pidc_pdf.groupby('A_pdbx_PDB_model_num'):
            pidc_coords = pdbid_chain_model[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
            if pidc_name == '1MSH_A' and model_num == 30:
                continue
            elif pidc_name == '2JSC_A' and model_num == 1:
                continue
            elif (pidc_name == '6UJV_A' or pidc_name == '6UJV_B' or
                  pidc_name == '6UJV_C' or pidc_name == '7CLV_A' or
                  pidc_name == '7CLV_B' or pidc_name == '8J4I_A' or
                  pidc_name == '8J4I_B' or pidc_name == '8J4I_C' or
                  pidc_name == '8J4I_D' or pidc_name == '8J4I_E' or
                  pidc_name == '8J4I_F'):
                continue
            else:
                rmsd = calculate_rmsd(coords1=mean_coords_pidc, coords2=pidc_coords)
                rmsds.append(round(rmsd, 4))
                pidc_model = f'{pidc_name}_{model_num}' if model_num >= 10 else f'{pidc_name}_0{model_num}'
                rmsd_per_model['pidc_model'].append(pidc_model)

        if len(rmsds) == 0:
            print(f'rmsds for {pidc_name} is empty: {rmsds}. Assigning np.nan.')
            min_rmsd, max_rmsd = np.nan, np.nan
            mean_rmsd, stdev_rmsd = np.nan, np.nan
        else:
            min_rmsd, max_rmsd = min(rmsds), max(rmsds)
            mean_rmsd, stdev_rmsd = statistics.mean(rmsds), statistics.stdev(rmsds)

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

        rmsd_pdf['pidc'] = [pidc_name] + ['' for _ in rmsds][: -1]
        rmsd_pdf['min_RMSD'] = [round(min_rmsd, 4)] + empty_list[: -1]
        rmsd_pdf['max_RMSD'] = [round(max_rmsd, 4)] + empty_list[: -1]
        rmsd_pdf['mean_RMSD'] = [round(mean_rmsd, 4)] + empty_list[: -1]
        rmsd_pdf['stdev_RMSD'] = [round(stdev_rmsd, 4)] + empty_list[: -1]

        rp_pidc_dir = os.path.join(_rp_rmsd_dir('multimod_2713_hetallchains_hom1chain'))
        os.makedirs(rp_pidc_dir, exist_ok=True)
        rmsd_pdf.to_csv(os.path.join(rp_pidc_dir, f'RMSD_{pidc_name}.csv'), index=False)

    print(f'Completed in {round((time() - start))} seconds.')  # 12 seconds


# Note: I don't superimpose coords of different models onto the randomly chosen ref model prior to computing the mean
# & stddev. It might be beneficial to do so, but for now I've opted not to.
def mean_stdev_struct(rp_pidchains_ssv: str, rp_dst_dir: str):
    pidchain = os.path.basename(rp_pidchains_ssv).removesuffix('.ssv')
    print(pidchain)
    pdf = pd.read_csv(rp_pidchains_ssv, sep=' ')
    # pdf = pdf.sort_values(by=["A_pdbx_PDB_model_num", "A_id"])
    model_numbers = pdf['A_pdbx_PDB_model_num'].unique()
    coord_list = []
    ref_A_id, ref_S_mon_id, ref_S_seq_id = None, None, None  # _atom_site.id is continous through all models (and chains).

    if pidchain == '1MSH_A':
        model_numbers = model_numbers[:-1]
        print("Excluding model 30 of 1MSH_A for now because missing 3 residues at C-term.")
    if pidchain == '2JSC_A':
        model_numbers = model_numbers[1:]
        print("Excluding model 1 of 2JSC_A for now because it lacks the first residue.")

    for i, model_num in enumerate(model_numbers):
        model_df = pdf[pdf['A_pdbx_PDB_model_num'] == model_num].copy()
        S_mon_id_values, S_seq_id_values = model_df['S_mon_id'].values, model_df['S_seq_id'].values

        if i == 0: # first model's labels for validation (asserts)
            ref_A_id, ref_S_mon_id, ref_S_seq_id = model_df['A_id'].values, S_mon_id_values, S_seq_id_values
        else:
            assert np.array_equal(S_mon_id_values, ref_S_mon_id), f'S_mon_id mismatch in model {model_num}'
            assert np.array_equal(S_seq_id_values, ref_S_seq_id), f'S_seq_id mismatch in model {model_num}'

        coords = model_df[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
        coord_list.append(coords)

    coord_array = np.stack(coord_list)  # Stack into 3D array: (n_models, n_atoms, 3)
    mean_coords, std_coords = coord_array.mean(axis=0), coord_array.std(axis=0, ddof=1)

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

    rp_rmsd_dst_dir = os.path.join(rp_dst_dir, 'mean_coords')
    os.makedirs(rp_rmsd_dst_dir, exist_ok=True)
    rp_rmsd_csv = os.path.join(rp_rmsd_dst_dir, f'{pidchain}.csv')
    if mean_df is not None:
        mean_df.to_csv(rp_rmsd_csv, index=False)
    return mean_df

if __name__ == '__main__':

    _calc_rmsds_and_stats()


    # rp_pidc_ssvs = sorted(glob.glob(os.path.join(_rp_parsed_cifs_dir('multimod_2713_hetallchains_hom1chain'), '*.ssv')))
    # rp_dst_dir = os.path.join(_rp_rmsd_dir('multimod_2713_hetallchains_hom1chain'))


    # for rp_pidc_ssv in rp_pidc_ssvs:
    #     mean_stdev_struct(rp_pidc_ssv, rp_dst_dir)
    #
    # rp_pidchains_ssvs = sorted(glob.glob(os.path.join('..', 'data', 'NMR', 'parsed_cifs',
    #                                            'multimod_2713_hetallchains_hom1chain', '*.ssv')))

    # rp_parsed_cifs_ssvs = sorted(glob.glob(os.path.join('..', 'data', 'NMR', 'parsed_cifs',
    #                                                     'multimod_2713_hetallchains_hom1chain', '*.ssv')))
    #
    # rp_mean_coords_csvs = sorted(glob.glob(os.path.join('..', 'data', 'NMR', 'RMSD',
    #                                                     'multimod_2713_hetallchains_hom1chain',
    #                                                     'mean_coords', '*.csv')))
    #
    # for rp_parsed_cifs_ssv_, rp_mean_coords_csv_ in zip(rp_parsed_cifs_ssvs, rp_mean_coords_csvs):
    #     rmsds_, model_nums_, pidch_list4table_ = calc_rmsds_of_models(rp_parsed_cifs_ssv=rp_parsed_cifs_ssv_,
    #                                                                   rp_mean_coords_csv=rp_mean_coords_csv_)


    # _meric = 'heteromeric'
    # with open(f'../data/NMR/multimodel_lists/{_meric[:3]}_multimod_2104_PidChains.txt', 'r') as f:
    #     pdbid_chains = f.readlines()
    # rmsd_mat, n_models = compute_rmsd_matrix(pdbid_chain='1A0N_A', het_hom=_meric)

    # # linkage_matrix = build_dendrogram(rmsd_mat, n_models)
    # # build_heatmap(rmsd_mat, linkage_matrix)
    # # build_contour_map(rmsd_mat)

    # cluster_save_to_json(rmsd_mat, pdbid_chain)
    #
    # # clusters = cluster_models(rmsd_mat, threshold=2.0)
    # # json_dict = clusters_to_json_dict(clusters)
    # # mm_ensbls = '../data/NMR/multimodel_ensembles'
    # # os.makedirs(mm_ensbls, exist_ok=True)
    # # output_file = f'{mm_ensbls}/{pdbid_chain}_ens.json'
    # # write_clusters_json(json_dict, output_file)
    # # print(f'Saved ensembles to {output_file}')

    # PATHS OF PDBid_CHAINS:
    # rp_mmseqs_results_dir = mmseqs2.rp_mmseqs_results_dir(het_hom='hethom_combined')    # TODO update
    # rp_homologues_30_20_90 = os.path.join(rp_mmseqs_results_dir, 'homologues_30_20_90.csv')
    # rp_non_homologues_30_20_90 = os.path.join(rp_mmseqs_results_dir, 'non_homologues_30_20_90.csv')

    # PidChains:
    # homologues_30_20_90_pdf = pd.read_csv(rp_homologues_30_20_90)
    # non_homologues_30_20_90_pdf = pd.read_csv(rp_non_homologues_30_20_90)

    # homologues_query = homologues_30_20_90_pdf['query'].tolist()
    # homologues_target = homologues_30_20_90_pdf['target'].tolist()
    # homologues_all = list(set(homologues_query + homologues_target))
    # homologues_all.sort()
    # homologues_pairs_row_list = list(zip(homologues_30_20_90_pdf['query'], homologues_30_20_90_pdf['target']))

    # non_homologues_query = non_homologues_30_20_90_pdf['query'].tolist()
    # non_homologues_target = non_homologues_30_20_90_pdf['target'].tolist()
    # non_homologues_all = list(set(non_homologues_query + non_homologues_target))
    # non_homologues_all.sort()
    #
    # all_pidChains = list(set(homologues_all + non_homologues_all))
    # all_pidChains.sort()  # NOTE: AFTER NOTICING THE TOTAL NUMBER OF PDBID_CHAINS IS LOWER THAN IT SHOULD BE,
    # # I FOUND OUT THE REASON IS THAT MMSEQS2 JUST DOESN'T BOTHER RETURNING ANY RESULTS FOR THOSE BELOW
    # # SOME INTERNAL THRESHOLD. SO I'M INSTEAD JUST USING THE FULL LISTS FROM THE MULTIMODEL PIDCHAINS LST FILES.
    # print(f'len(all_pidChains)={len(all_pidChains)}')

    # with open(f'../data/NMR/multimodel_lists/het_multimod_2104_PidChains.txt', 'r') as f:
    #     het_pidChains = f.read().split()
    # with open(f'../data/NMR/multimodel_lists/hom_multimod_1421_PidChains.txt', 'r') as f:
    #     hom_pidChains = f.read().split()
    #
    # all_pidChains_from_lists = list(set(het_pidChains + hom_pidChains))
    # all_pidChains_from_lists.sort()
    # print(f'len(all_pidChains_from_lists)={len(all_pidChains_from_lists)}')

    # rp_parsed_cif_dir = os.path.join('..', 'data', 'NMR', 'parsed_cifs', 'hethom_combined')
    #


    # rmsd_results, pidChain1_list, pidChain2_list = [], [], []
    # pidChain_rmsds = []
    # rp_rmsd_mean_coords_dir = os.path.join('..', 'data', 'NMR', 'RMSD', 'hethom_combined', 'mean_coords')
    # for pidChain1_, pidChain2_ in homologues_pairs_row_list:
    #     rmsd_result = calc_rmsds_pidchain_pairs(pidChain1_, pidChain2_, rp_rmsd_mean_coords_dir)
    #     pidChain1_list.append(pidChain1_)
    #     pidChain2_list.append(pidChain2_)
    #
    #     pidChain_rmsds.append({
    #         'query': pidChain1_,
    #         'target': pidChain2_,
    #         'RMSD': rmsd_result,
    #     })
    #
    # pdf = pd.DataFrame(pidChain_rmsds)
    # pdf_sorted = pdf.sort_values(by=['RMSD'], ascending=[True])
    # rp_rmsd_dir = _rp_RMSD_dir(het_hom='hethom_combined')
    # os.makedirs(rp_rmsd_dir, exist_ok=True )
    # pdf.to_csv(os.path.join(rp_rmsd_dir, 'homologous_30_20_90.csv'), index=False)

