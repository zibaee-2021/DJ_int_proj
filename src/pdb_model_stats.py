import os, glob
from time import time
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
import RMSD
import mmseqs2


def _total_chain_count_and_year(het_hom: str, pdbid: str, use_mmcif: bool) -> tuple:
    parser = MMCIFParser(QUIET=True)
    if not use_mmcif:
        parser = PDBParser(QUIET=True)
    rp_raw_cif_dir = os.path.join('..', 'data', 'NMR', 'raw_cifs', het_hom)
    rp_cif = os.path.join(rp_raw_cif_dir, f'{pdbid}.cif')
    bio_struct = parser.get_structure('', rp_cif)
    total_chain_count = len(list(next(bio_struct.get_models()).get_chains()))
    year = bio_struct.header['deposition_date'][:4]
    return total_chain_count, year


def pdbid_dict_chain(rp_pdbid_chains: str) -> dict:
    with open(rp_pdbid_chains, 'r') as f:
        pdbid_chains = [line.rstrip('\n') for line in f]
    # 2104 heteromeric PDBid_chains
    pdbid_chains.sort()  # expected to already be sorted (i.e. before writing out), but to avoid relying on other code.
    pdbid_chains_dict = defaultdict(list)

    for pdbid_chain in pdbid_chains:
        pid, ch = pdbid_chain.split('_')
        pdbid_chains_dict[pid].append(ch)
    # 932 heteromeric PDBids
    return dict(pdbid_chains_dict)


def _read_multimodel_pdbid_chains(txt_f: str):
    rp_pdbid_chains = os.path.join('..', 'data', 'NMR', 'multimodel_lists', txt_f)
    pdbid_chains_dict = pdbid_dict_chain(rp_pdbid_chains)
    # 2104 heteromeric pdbid_chains --> 932 heteromeric pdbids.
    pdbids = list(pdbid_chains_dict.keys())
    return pdbids, pdbid_chains_dict


def _calc_rmsds(pdbid_chain: str, rp_token_cif_dir: str) -> tuple:
    # rmsd_mat, n_models = RMSD.compute_rmsd_matrix(pdbid_chain, het_hom)
    # clusters = RMSD.cluster_models(rmsd_mat, threshold=2.0)
    ref_structure = RMSD.mean_stddev_struct(pdbid_chain)
    np_ref_model_coords = ref_structure[['mean_x', 'mean_y', 'mean_z']].values
    rp_pdbid_chain_ssv = os.path.join(rp_token_cif_dir, f'{pdbid_chain}.ssv')
    pdbid_chain_pdf = pd.read_csv(rp_pdbid_chain_ssv, sep=' ')
    rmsds, model_nums = [], []
    pdbidchain_list = [pdbid_chain]
    for model_num, pdbid_chain_model in pdbid_chain_pdf.groupby('A_pdbx_PDB_model_num'):
        np_pdbidchain_coords = pdbid_chain_model[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
        rmsd_pdbid_chain = RMSD.rmsd_reference(np_model1_coords=np_ref_model_coords,
                                               np_model2_coords=np_pdbidchain_coords)
        model_nums.append(model_num)
        rmsds.append(rmsd_pdbid_chain)
        pdbidchain_list.append('')  # so that the resulting stats table has the pdbid just on the first row only.
    pdbidchain_list = pdbidchain_list[:-1]
    rmsds = np.array(rmsds, dtype=np.float16)
    return rmsds, model_nums, pdbidchain_list


def _calc_identity_for_stats(het_hom: str, pdbid: str, filt_pdf) -> str:
    if not filt_pdf.empty:
        print(f'There is some heterogeneity of sequence between some chains. '
              f'See data/NMR/mmseqs/{het_hom}/results/{pdbid}.csv')
        identity = '(chains not identical, see mmseqs pdf)'
    else:
        print('All chains have high identity, if not 100%')
        identity = '0'
    return identity


def _calc_mmseqs2(het_hom: str, pid: str, chains: list):
    # rp_fasta_file = mmseqs2.write_fasta(het_hom, pid, chains)
    # rp_fasta_file = mmseqs2.combine_all_fastas_in_1_file(het_hom)
    # result_pdf = mmseqs2.run_mmseqs2_all_vs_all(rp_fasta_file)

    # rp_dst_combo_fastas_f = os.path.join(mmseqs2._rp_mmseqs_comb_fastas_dir(het_hom), 'het_all_932_PDBids.fasta')
    # if het_hom == 'homomeric':
    #     rp_dst_combo_fastas_f = os.path.join(mmseqs2._rp_mmseqs_comb_fastas_dir(het_hom), 'hom_all_609_PDBids.fasta')

    # result_pdf = mmseqs2.run_easysearch_mmseqs2_allvsall(rp_dst_combo_fastas_f)
    # result_pdf = mmseqs2.run_easysearch_mmseqs2_allvsall(rp_fasta_file)
    # filtered_pdf = mmseqs2.filter_and_write_results(het_hom, pid, result_pdf)

    result_pdf = mmseqs2.run_mmseqs2_across_all_pdbids_combined()
    filtered_pdf = mmseqs2.filter_and_write_results(het_hom='', pdbid='all_1541_hethom', pdf=result_pdf)

    return filtered_pdf


def generate_stats_het(het_hom: str, pdbid_chains_txt: str, use_mmcif=True):
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

    4. Calculate sequence identity between every chain for each PDBid from list, by MMseqs2.

    5. Calculate RMSD for each model against the average of those models.
    """
    start = time()
    stats = []
    rp_token_cif_dir = os.path.join('..', 'data', 'NMR', 'tokenised_cifs', het_hom)
    pdbids, pdbid_chains_dict = _read_multimodel_pdbid_chains(txt_f=pdbid_chains_txt)

    for pid, chains in pdbid_chains_dict.items():
        total_chain_count, year = _total_chain_count_and_year(het_hom, pid, use_mmcif)
        mmseqs_pdf = _calc_mmseqs2(het_hom, pid, chains)

        for chain in chains:
            identity = _calc_identity_for_stats(het_hom, pid, mmseqs_pdf)
            rmsds, model_nums, pdbidchain_list = _calc_rmsds(f'{pid}_{chain}', rp_token_cif_dir)
            max_rmsd, min_rmsd = np.max(rmsds), np.min(rmsds)
            rmsd_pdf = pd.DataFrame({'pdbid_chain': pdbidchain_list, 'model_num':model_nums, 'rmsd': rmsds})
            rp_rmsd_dir = os.path.join('..', 'data', 'NMR', 'RMSD', het_hom)
            os.makedirs(rp_rmsd_dir, exist_ok=True)
            rp_rmsd_csv = os.path.join(rp_rmsd_dir, f'{pid}_{chain}.csv')
            rmsd_pdf.to_csv(rp_rmsd_csv, index=False)
            rp_tok_cif = os.path.join(rp_token_cif_dir, f'{pid}_{chain}.ssv')
            pdbid_chain_pdf = pd.read_csv(rp_tok_cif, sep=' ')
            model_count = len(pdbid_chain_pdf['A_pdbx_PDB_model_num'].unique())
            ca_count = pdbid_chain_pdf.shape[0] / model_count

            stats.append({
                'PDBid': pid,
                'het_hom': het_hom,
                'year': year,
                '#model': model_count,
                '#chains': total_chain_count,
                '#protChains': len(pdbid_chains_dict[pid]),
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
    print(f'Completed 2104 CIFs in {round((time() - start)/60)} mins')
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


if __name__ == '__main__':

    # _calc_mmseqs2(het_hom='heteromeric', pid='0_all_het', chains=[])
    # _calc_mmseqs2(het_hom='homomeric', pid='0_all_hom', chains=[])
    _calc_mmseqs2(het_hom='', pid='', chains=[])

    # _calc_mmseqs2(het_hom='heteromeric', pdbid='1AI0',
    #               chains=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'])
    # rpath_pdbid_chains = (f'../data/NMR/multimodel_lists/het_multimod_2104_pdbid_chains.txt')
    # stats_pdf = generate_stats_het(het_hom='heteromeric', pdbid_chains_txt='het_multimod_2104_pdbid_chains.txt')
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


    # Note:
    # cif_dict = MMCIF2Dict.MMCIF2Dict(rpath_cif)
    # A few of these return '?', so using 'get_structure()' instead,
    # which broadly same but sometimes differs by as much as 3 years.
    # year = cif_dict.get('_citation.year', ['NA'])[0]
    # if year != 'NA' and year != '?':
    #     print(f'year is given as {year}. Replacing with 0')
    #     year = int(year)
    # else:
    #     year = 0

    # from Bio.PDB.Polypeptide import is_aa
    # first_model = next(structure.get_models())
    # for chain in first_model:
    #     ca_count = sum(1 for atom in chain.get_atoms() if atom.get_name() == 'CA')
    # num_models = len(list(structure))