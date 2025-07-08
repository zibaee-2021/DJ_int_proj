import os, glob, json
from time import time
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

