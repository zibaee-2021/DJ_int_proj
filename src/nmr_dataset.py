import os, glob, json
from time import time
import pandas as pd


from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from src.utils import api_callr as api
from src.utils import cif_parsr as cp

# start = time()
# relpath_raw_cifs = os.path.join('..', 'data', 'NMR', 'raw_cifs')
# _meric = 'homomeric'
# relpath_raw_cifs_meric = os.path.join(relpath_raw_cifs, _meric)
# os.makedirs(relpath_raw_cifs_meric, exist_ok=True)
# sol_nmr_homo_686_pdbids = api.call_rcsb_for_pdbids_of_solution_nmr_homomeric(number=686)
#
# for pdbid in sol_nmr_homo_686_pdbids:
#     response = api.call_rcsb_for_cif(pdbid)
#     with open(os.path.join(relpath_raw_cifs_meric, f'{pdbid}.cif'), 'w') as cif_file:
#         cif_file.write(response.text)
#
# _meric = 'heteromeric'
# relpath_raw_cifs_meric = os.path.join(relpath_raw_cifs, _meric)
# os.makedirs(relpath_raw_cifs_meric, exist_ok=True)
# sol_nmr_hetero_1038_pdbids = api.call_rcsb_for_pdbids_of_solution_nmr_heteromeric(number=1038)
#
# for pdbid in sol_nmr_hetero_1038_pdbids:
#     response = api.call_rcsb_for_cif(pdbid)
#     with open(os.path.join(relpath_raw_cifs_meric, f'{pdbid}.cif'), 'w') as cif_file:
#         cif_file.write(response.text)
# print(f'Completed {len(sol_nmr_homo_686_pdbids)} homomeric PDBs and {len(sol_nmr_hetero_1038_pdbids)} heteromeric PDBs in {round((time() - start) / 60)} minutes')


# # 2. Parse 686 homomeric proteins, write to ssvs:
# # 3. Parse 1038 heteromeric proteins, write to ssvs:
# import os, glob
# from time import time
# from Bio.PDB.MMCIF2Dict import MMCIF2Dict
# from src.utils import cif_parsr as cp
# start = time()
# # _meric = 'homomeric' # takes about 80 minutes
# _meric = 'heteromeric'  # takes about 140 minutes
# relpath_raw_cifs_meric = os.path.join('..', 'data', 'NMR', 'raw_cifs', _meric)
# relpath_token_cifs = os.path.join('..', 'data', 'NMR', 'tokenised_cifs', _meric)
# os.makedirs(relpath_token_cifs, exist_ok=True)
#
# relpath_cifs = glob.glob(os.path.join(relpath_raw_cifs_meric, f'*.cif'))
# no_CA_pdbids = list()
#
# for relpath_cif in relpath_cifs:
#     cif_dict = MMCIF2Dict(relpath_cif)
#     pdbid = os.path.basename(relpath_cif).removesuffix('.cif')
#     cif_pdfs_per_chain = cp.parse_cif(pdb_id=pdbid, mmcif_dict=cif_dict)  # same function from MSc project, but alpha-carbons only.
#     # with open(os.path.join('..', 'data', 'enumeration', 'residues.json'), 'r') as json_f:
#     #     residues_enumerated = json.load(json_f)
#
#     # for pdf_chain in cif_pdfs_per_chain:
#     #     chain = pdf_chain['S_asym_id'].iloc[0]
#         # pdf_chain = pdf_chain.copy()
#         # pdf_chain.loc[:, 'aa_label_num'] = pdf_chain['S_mon_id'].map(residues_enumerated).astype('Int64')
#         # pdf_chain = pdf_chain[['A_pdbx_PDB_model_num', 'S_seq_id', 'S_mon_id', 'aa_label_num',
#         #                        'A_id', 'A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']]
#         # pdf_chain.to_csv(path_or_buf=os.path.join(relpath_token_cifs, f'{pdbid}_{chain}.ssv'), sep=' ', index=False)
#
# print(f'Completed {len(relpath_cifs)} {_meric} PDBs in {round((time() - start) / 60)} minutes')

# # 2. Parse 686 homomeric proteins, write to ssvs:
# # 3. Parse 1038 heteromeric proteins, write to ssvs:
# import os, glob
# from time import time
# from Bio.PDB.MMCIF2Dict import MMCIF2Dict
# from src.utils import cif_parsr as cp
# start = time()
# # _meric = 'homomeric' # takes about 80 minutes
# _meric = 'heteromeric'  # takes about 140 minutes
# relpath_raw_cifs_meric = os.path.join('..', 'data', 'NMR', 'raw_cifs', _meric)
# relpath_token_cifs = os.path.join('..', 'data', 'NMR', 'tokenised_cifs', _meric)
# os.makedirs(relpath_token_cifs, exist_ok=True)
#
# relpath_cifs = glob.glob(os.path.join(relpath_raw_cifs_meric, f'*.cif'))
# no_CA_pdbids = list()
#
# for relpath_cif in relpath_cifs:
#     cif_dict = MMCIF2Dict(relpath_cif)
#     pdbid = os.path.basename(relpath_cif).removesuffix('.cif')
#     cif_pdfs_per_chain = cp.parse_cif(pdb_id=pdbid, mmcif_dict=cif_dict)  # same function from MSc project, but alpha-carbons only.
#     # with open(os.path.join('..', 'data', 'enumeration', 'residues.json'), 'r') as json_f:
#     #     residues_enumerated = json.load(json_f)
#
#     # for pdf_chain in cif_pdfs_per_chain:
#     #     chain = pdf_chain['S_asym_id'].iloc[0]
#         # pdf_chain = pdf_chain.copy()
#         # pdf_chain.loc[:, 'aa_label_num'] = pdf_chain['S_mon_id'].map(residues_enumerated).astype('Int64')
#         # pdf_chain = pdf_chain[['A_pdbx_PDB_model_num', 'S_seq_id', 'S_mon_id', 'aa_label_num',
#         #                        'A_id', 'A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']]
#         # pdf_chain.to_csv(path_or_buf=os.path.join(relpath_token_cifs, f'{pdbid}_{chain}.ssv'), sep=' ', index=False)
#
# print(f'Completed {len(relpath_cifs)} {_meric} PDBs in {round((time() - start) / 60)} minutes')


# # 2. Parse 686 homomeric proteins, write to ssvs:
# # 3. Parse 1038 heteromeric proteins, write to ssvs:
# import os, glob
# from time import time
# from Bio.PDB.MMCIF2Dict import MMCIF2Dict
# from src.utils import cif_parsr as cp
# start = time()
# # _meric = 'homomeric' # takes about 80 minutes
# _meric = 'heteromeric'  # takes about 140 minutes
# relpath_raw_cifs_meric = os.path.join('..', 'data', 'NMR', 'raw_cifs', _meric)
# relpath_token_cifs = os.path.join('..', 'data', 'NMR', 'tokenised_cifs', _meric)
# os.makedirs(relpath_token_cifs, exist_ok=True)
#
# relpath_cifs = glob.glob(os.path.join(relpath_raw_cifs_meric, f'*.cif'))
# no_CA_pdbids = list()
#
# for relpath_cif in relpath_cifs:
#     cif_dict = MMCIF2Dict(relpath_cif)
#     pdbid = os.path.basename(relpath_cif).removesuffix('.cif')
#     cif_pdfs_per_chain = cp.parse_cif(pdb_id=pdbid, mmcif_dict=cif_dict)  # same function from MSc project, but alpha-carbons only.
#     # with open(os.path.join('..', 'data', 'enumeration', 'residues.json'), 'r') as json_f:
#     #     residues_enumerated = json.load(json_f)
#
#     # for pdf_chain in cif_pdfs_per_chain:
#     #     chain = pdf_chain['S_asym_id'].iloc[0]
#         # pdf_chain = pdf_chain.copy()
#         # pdf_chain.loc[:, 'aa_label_num'] = pdf_chain['S_mon_id'].map(residues_enumerated).astype('Int64')
#         # pdf_chain = pdf_chain[['A_pdbx_PDB_model_num', 'S_seq_id', 'S_mon_id', 'aa_label_num',
#         #                        'A_id', 'A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']]
#         # pdf_chain.to_csv(path_or_buf=os.path.join(relpath_token_cifs, f'{pdbid}_{chain}.ssv'), sep=' ', index=False)
#
# print(f'Completed {len(relpath_cifs)} {_meric} PDBs in {round((time() - start) / 60)} minutes')


# THIS (AND NEXT) CELL USE THE PDBS FROM THE TOKENISED DIR, HENCE ALREADY FILTERED OUT THOSE WITH NO CA ATOMS:
if __name__ == "__main__":
    start = time()
    # _meric = 'homomeric'  # takes about 8 seconds
    _meric = 'heteromeric'  # takes about 16 seconds
    relpath_token_cifs = os.path.join('..', 'data', 'NMR', 'tokenised_cifs', _meric)

    ssv_files = glob.glob(os.path.join(relpath_token_cifs, '*.ssv'))

    multimodel_pdbids = []
    single_model_pdbids = []

    for file_path in ssv_files:
        try:
            df = pd.read_csv(file_path, sep=' ', dtype=str)
            unique_models = df['A_pdbx_PDB_model_num'].unique()
            if len(unique_models) > 1:
                multimodel_pdbids.append(os.path.basename(file_path))
            else:
                single_model_pdbids.append(os.path.basename(file_path))
        except Exception as e:
            print(f'Error reading {file_path}: {e}')

        single_model_pdbids.sort()
        multimodel_pdbids.sort()

    list_dir = os.path.join('..', 'data', 'NMR', 'multimodel_PDBids')
    os.makedirs(list_dir, exist_ok=True)

    singlemodel_txt = os.path.join(list_dir, f'{_meric[:3]}_singlemod_{len(single_model_pdbids)}_pdbid_chains.txt')
    with open(singlemodel_txt, 'w') as f:
        headr = f'The following {len(single_model_pdbids)} {_meric} PDBid_chains had only 1 model (solution NMR) in the RCSB dataset:\n'
        singlemods = list()
        for i, fname in enumerate(single_model_pdbids):
            j = i + 1
            newline = '\n' if j % 16 == 0 else ''
            singlemods.append(f'{fname.removesuffix('.ssv')} {newline}')
        singlemods = headr + ''.join(singlemods)
        f.write(singlemods)

    multimodel_txt = os.path.join(list_dir, f'{_meric[:3]}_multimod_{len(multimodel_pdbids)}_pdbid_chains.txt')
    with open(multimodel_txt, 'w') as f:
        for fname in multimodel_pdbids:
            f.write(fname.removesuffix('.ssv') + '\n')

    print(f'\nSaved {len(multimodel_pdbids)} {_meric} pdbid_chains with >1 unique model to text file.')
    print(f'Completed in {round(time() - start)} seconds')
