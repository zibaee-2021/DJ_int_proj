import os, glob, json
from time import time
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from src.utils import cif_parsr as cp
from src.utils import api_callr as api


def write_struct_files_for_solution_NMR(_meric: str, cif_or_pdb: str) -> None:
    start = time()
    dst_dir = os.path.join('..', 'data', 'NMR', f'raw_{cif_or_pdb}s', _meric)
    os.makedirs(dst_dir, exist_ok=True)
    if _meric == 'homomeric':
        sol_nmr_pdbids = api.call_rcsb_for_pdbids_of_solution_nmr_homomeric(number=686)
    else:
        sol_nmr_pdbids = api.call_rcsb_for_pdbids_of_solution_nmr_heteromeric(number=1038)

    for pdbid in sol_nmr_pdbids:
        response = api.call_rcsb_for_cif_or_pdb(pdbid, cif_or_pdb=cif_or_pdb)
        with open(os.path.join(dst_dir, f'{pdbid}.{cif_or_pdb}'), 'w') as f:
            f.write(response.text)
    print(f'Completed {len(sol_nmr_pdbids)} {_meric} {cif_or_pdb}s in {round((time() - start) / 60)} minutes')


def parse_atomic_records_from_cifs(_meric: str, write_results=True):
    start = time()
    rp_raw_cifs_meric = os.path.join('..', 'data', 'NMR', 'raw_cifs', _meric)
    rp_token_cifs = os.path.join('..', 'data', 'NMR', 'parsed_cifs', _meric)
    os.makedirs(rp_token_cifs, exist_ok=True)

    rpath_cifs = glob.glob(os.path.join(rp_raw_cifs_meric, f'*.cif'))
    rpath_cifs.sort()

    for rp_cif in rpath_cifs:
        cif_dict = MMCIF2Dict(rp_cif)
        pdbid = os.path.basename(rp_cif).removesuffix('.cif')
        # pdbid = '2N2K'
        cif_pdfs_per_chain, empty_pdbidchains = cp.parse_cif(pdb_id=pdbid, mmcif_dict=cif_dict)  # same function from MSc project, but alpha-carbons only.
        empty_pdbidchains.sort()

        for pdbid_chain in empty_pdbidchains:
            with open(os.path.join('..', 'data', 'NMR', 'tokenised_cifs', f'noCA_pdbidchains_{_meric[:3]}.txt'), 'a') as f:
                f.write(pdbid_chain + '\n')

        with open(os.path.join('..', 'data', 'enumeration', 'residues.json'), 'r') as json_f:
            residues_enumerated = json.load(json_f)
        for pdf_chain in cif_pdfs_per_chain:
            try:
                chain = pdf_chain['S_asym_id'].iloc[0]
            except:
                print(f"pdf_chain['S_asym_id'].iloc[0] on line 49 is failing for {pdbid}")
                print(f"pdf_chain={pdf_chain}")
                print('Leaving this one out.')
                continue
            pdf_chain = pdf_chain.copy()
            pdf_chain.loc[:, 'aa_label_num'] = pdf_chain['S_mon_id'].map(residues_enumerated).astype('Int64')
            pdf_chain = pdf_chain[['A_pdbx_PDB_model_num', 'S_seq_id', 'S_mon_id', 'aa_label_num',
                                   'A_id', 'A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']]
            if write_results:
                pdf_chain.to_csv(path_or_buf=os.path.join(rp_token_cifs, f'{pdbid}_{chain}.ssv'), sep=' ', index=False)

    empty_pdbidchains_all.sort()
    if write_results:
        for pdbid_chain in empty_pdbidchains_all:
            with open(os.path.join('..', 'data', 'NMR', 'parsed_cifs', f'{_meric[:3]}_noCA_PidChains.txt'), 'a') as f:
                f.write(pdbid_chain + '\n')

    print(f'Completed {len(rpath_cifs)} {_meric} CIFs in {round((time() - start) / 60)} minutes')


def generate_pdb_lists_from_parsed_ssvs(_meric: str):
    start = time()
    rp_token_cifs = os.path.join('..', 'data', 'NMR', 'parsed_cifs', _meric)
    ssv_files = glob.glob(os.path.join(rp_token_cifs, '*.ssv'))

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

    list_dir = os.path.join('..', 'data', 'NMR', 'multimodel_lists')
    os.makedirs(list_dir, exist_ok=True)

    singlemodel_txt = os.path.join(list_dir, f'{_meric[:3]}_singlemod_{len(single_model_pdbids)}_PidChains.txt')
    with open(singlemodel_txt, 'w') as f:
        headr = f'The following {len(single_model_pdbids)} {_meric} PDBid_chains had only 1 model (solution NMR) in the RCSB dataset:\n'
        singlemods = list()
        for i, fname in enumerate(single_model_pdbids):
            j = i + 1
            newline = '\n' if j % 16 == 0 else ''
            singlemods.append(f'{fname.removesuffix('.ssv')} {newline}')
        singlemods = headr + ''.join(singlemods)
        f.write(singlemods)

    multimodel_txt = os.path.join(list_dir, f'{_meric[:3]}_multimod_{len(multimodel_pdbids)}_PidChains.txt')
    with open(multimodel_txt, 'w') as f:
        for fname in multimodel_pdbids:
            f.write(fname.removesuffix('.ssv') + '\n')

    print(f'\nSaved {len(multimodel_pdbids)} {_meric} pdbid_chains with >1 unique model to text file.')
    print(f'Completed in {round(time() - start)} seconds')


def write_pdb_or_cif(pdbid: str, cif_or_pdb: str, dst_dir: str):
    import os
    from src.utils import api_callr as api
    response = api.call_rcsb_for_cif_or_pdb(pdb_id=pdbid, cif_or_pdb=cif_or_pdb)
    with open(os.path.join(dst_dir, f'{pdbid}.{cif_or_pdb}'), 'w') as f:
        f.write(response.text)



if __name__ == "__main__":
    _meric = 'homomeric'
    # write_struct_files_for_solution_NMR(_meric=_meric, cif_or_pdb='pdb')  # 16 mins (new Mac) for 686-4 = 682 PDBs. (21 mins Rocky)
    # write_struct_files_for_solution_NMR(_meric=_meric, cif_or_pdb='cif')  # 21 mins (new Mac) for 686 mmCIFs. (21 mins Rocky)
    parse_atomic_records_from_cifs(_meric=_meric, write_results=True)  # 6 mins to parse (new Mac). (13 mins Rocky)
    # generate_pdb_lists_from_parsed_ssvs(_meric=_meric)  # 2 secs to generate PDBid_chain lists (multi & single model) (4 secs Rocky)

    _meric = 'heteromeric'
    write_struct_files_for_solution_NMR(_meric=_meric, cif_or_pdb='pdb')  # 31 mins (new Mac) for 1038 PDBs.  (31 mins Rocky)
    write_struct_files_for_solution_NMR(_meric=_meric, cif_or_pdb='cif')  # 32 mins (new Mac) for 1038 mmCIFs. (32 mins Rocky)
    parse_atomic_records_from_cifs(_meric=_meric)  # 8 mins to parse (new Mac). (17 mins Rocky)
    generate_pdb_lists_from_parsed_ssvs(_meric=_meric)  # 4 seconds to generate PDBid_chain lists (multi & single model)

# Notice: 4 legacy PDB could not be read (7ZE0, 9D9A, 9D9B, 9D9C):
# Failed to retrieve data from API: 404 Client Error: Not Found for url: https://files.rcsb.org/download/7ZE0.pdb
# Failed to retrieve data from API: 404 Client Error: Not Found for url: https://files.rcsb.org/download/9D9A.pdb
# Failed to retrieve data from API: 404 Client Error: Not Found for url: https://files.rcsb.org/download/9D9B.pdb
# Failed to retrieve data from API: 404 Client Error: Not Found for url: https://files.rcsb.org/download/9D9C.pdb

    # cif_or_pdb = 'cif'
    # dst_dir = os.path.join('..', 'data', 'NMR', f'raw_{cif_or_pdb}s', 'homomeric')
    # os.makedirs(dst_dir, exist_ok=True)
    # write_pdb_or_cif('7ZE0', cif_or_pdb, dst_dir)
    # write_pdb_or_cif('9D9A', cif_or_pdb, dst_dir)
    # write_pdb_or_cif('9D9B', cif_or_pdb, dst_dir)
    # write_pdb_or_cif('9D9C', cif_or_pdb, dst_dir)

# Possible errors during parsing:
# parsing 6F3K
# pdf_chain['S_asym_id'].iloc[0] on line 49 is failing for 6F3K
# pdf_chain=Empty DataFrame
# Columns: [A_pdbx_PDB_model_num, S_asym_id, S_seq_id, S_mon_id, A_id, A_label_atom_id, A_Cartn_x, A_Cartn_y, A_Cartn_z]
# Index: []
# Leaving this one out.

# parsing 2N2K
# pdf_chain['S_asym_id'].iloc[0] on line 49 is failing for 2N2K
# pdf_chain=Empty DataFrame
# Columns: [A_pdbx_PDB_model_num, S_asym_id, S_seq_id, S_mon_id, A_id, A_label_atom_id, A_Cartn_x, A_Cartn_y, A_Cartn_z]
# Index: []
# Leaving this one out.

# The following PDBs included chains with non-natural amino acids, so `_atom_site` which will be filtered to
# remove HETATM will have fewer chains than the corresponding `_pdbx_poly_seq_scheme`:
# parsing 2KMJ. 1 chains in _atom_site, 3 chains in _pdbx_poly_seq_scheme.
# parsing 1S4A. 0 chains in _atom_site, 2 chains in _pdbx_poly_seq_scheme.
# parsing 1S1O. 0 chains in _atom_site, 2 chains in _pdbx_poly_seq_scheme.
#
# parsing 7B3K. 1 chains in _atom_site, 2 chains in _pdbx_poly_seq_scheme.
# parsing 7B3J. 1 chains in _atom_site, 2 chains in _pdbx_poly_seq_scheme.