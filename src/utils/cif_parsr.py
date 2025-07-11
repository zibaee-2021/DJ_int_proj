"""
'A_' = `_atom_site`.
'S_' = `_pdbx_poly_seq_scheme`.

A_group_PDB       # 'ATOM' or 'HETATM'    - FILTER ON THIS.                                           THEN REMOVE
S_pdb_seq_num,    # RESIDUE POSITION      - (NOT USED)                                                JUST REMOVE
A_label_seq_id,   # RESIDUE POSITION      - USED TO JOIN WITH S_seq_id.                               THEN REMOVE
A_label_comp_id,  # RESIDUE (3-LETTER)    - USED TO SANITY-CHECK WITH S_mon_id.                       THEN REMOVE
A_label_asym_id,  # CHAIN                 - JOIN ON THIS, SORT ON THIS.                               THEN REMOVE
A_occupancy       # OCCUPANCY             - FILTER ON THIS.                                           THEN REMOVE

S_asym_id,        # CHAIN                 - JOIN ON THIS.                                             THEN KEEP
S_seq_id,         # RESIDUE POSITION      - USED TO JOIN WITH A_label_seq_id. SORT ON THIS.           THEN KEEP
S_mon_id,         # RESIDUE (3-LETTER)    - USE FOR SANITY-CHECK AGAINST A_label_comp_id.             THEN KEEP
A_id,             # ATOM POSITION         - SORT ON THIS.                                             THEN KEEP
A_label_atom_id,  # ATOM                  -                                                           JUST KEEP
A_Cartn_x,        # COORDINATES           - ATOM X-COORDINATES                                        JUST KEEP
A_Cartn_y,        # COORDINATES           - ATOM Y-COORDINATES                                        JUST KEEP
A_Cartn_z,        # COORDINATES           - ATOM Z-COORDINATES                                        JUST KEEP
"""

import os
from typing import List, Tuple
import numpy as np
import pandas as pd



def _impute_missing_coords(pdf_to_impute, value_to_impute_with=0):
    missing_count = pdf_to_impute['A_Cartn_x'].isna().sum()
    print(f"BEFORE imputing, {missing_count} rows with missing values in column 'A_Cartn_x'")
    pdf_to_impute[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_x']] = (pdf_to_impute[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']]
                                                              .fillna(value_to_impute_with, inplace=False))

    missing_count = pdf_to_impute['A_Cartn_x'].isna().sum()
    assert missing_count == 0, (f'AFTER imputing, there should be no rows with missing values, '
                                f"but {missing_count} rows in column 'A_Cartn_x' have NANs. "
                                f'Therefore something has gone wrong!')
    return pdf_to_impute


def _process_missing_data(pdf_with_missing_data: pd.DataFrame, impute=False) -> pd.DataFrame:
    if impute:
        result_pdf = _impute_missing_coords(pdf_with_missing_data)
    else:
        result_pdf = pdf_with_missing_data.dropna(how='any', axis=0, inplace=False, subset=['A_Cartn_x'])
    return result_pdf


def _replace_low_occupancy_coords_with_nans(pdf: pd.DataFrame) -> pd.DataFrame:
    pdf['A_Cartn_x'] = np.where(pdf['A_occupancy'] <= 0.5, np.nan, pdf['A_Cartn_x'])
    pdf['A_Cartn_y'] = np.where(pdf['A_occupancy'] <= 0.5, np.nan, pdf['A_Cartn_y'])
    pdf['A_Cartn_z'] = np.where(pdf['A_occupancy'] <= 0.5, np.nan, pdf['A_Cartn_z'])
    return pdf


def _sort_by_chain_residues_atoms(pdf: pd.DataFrame) -> pd.DataFrame:
    # AS `pdf` WILL ONLY BE PASSED HERE FOR ONE CHAIN AT A TIME,
    # I ONLY NEED TO SORT ROWS BY MODEL NUMBER, RESIDUE SEQUENCE NUMBERING (SEQ ID) THEN ATOM SEQUENCE NUMBERING (A_ID):
    pdf.reset_index(drop=True, inplace=True)
    pdf = pdf.sort_values(['A_pdbx_PDB_model_num', 'S_seq_id', 'A_id'], ignore_index=True)
    return pdf


def _cast_objects_to_stringdtype(pdf: pd.DataFrame) -> pd.DataFrame:
    """
    pdf MUST BE FOR ONE POLYPEPTIDE CHAIN ONLY.
    """
    cols_to_cast = ['S_mon_id', 'A_label_comp_id', 'A_label_atom_id', 'A_label_asym_id' , 'S_asym_id']
    for col_to_cast in cols_to_cast:
        pdf[col_to_cast] = pdf[col_to_cast].astype('string')
    return pdf


def _cast_number_strings_to_numeric_types(pdf_merged: pd.DataFrame) -> pd.DataFrame:
    """
    pdf_merged MUST BE FOR ONE POLYPEPTIDE CHAIN ONLY.
    """
    # CAST STRINGS OF FLOATS TO NUMERIC:
    # for col in ['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z', 'A_occupancy', 'b_iso_or_equiv']:
    for col in ['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z', 'A_occupancy']:
        pdf_merged.loc[:, col] = pd.to_numeric(pdf_merged[col], errors='coerce')

    # CAST STRINGS OF INTS TO NUMERIC AND THEN TO INTEGERS:
    list_of_cols_to_cast = ['S_seq_id', 'A_label_seq_id', 'S_pdb_seq_num', 'A_id', 'A_pdbx_PDB_model_num']

    for col_to_cast in list_of_cols_to_cast:
        pdf_merged.loc[:, col_to_cast] = pd.to_numeric(pdf_merged[col_to_cast], errors='coerce')
        pdf_merged.loc[:, col_to_cast] = pdf_merged[col_to_cast].astype('Int64')
    return pdf_merged


def _rearrange_cols(pdf_merged: pd.DataFrame) -> pd.DataFrame:
    """
    pdf_merged MUST BE FOR ONE POLYPEPTIDE CHAIN ONLY.
    """
    return pdf_merged[[
        'A_pdbx_PDB_model_num',     # MODEL NUMBER
        'S_seq_id',                 # RESIDUE POSITION      - JOIN TO 'S_label_seq_id', SORT ON THIS, KEEP IN DF.
        'A_label_seq_id',           # RESIDUE POSITION      - JOIN TO 'S_seq_id', THEN REMOVE IT.
        'S_pdb_seq_num',            # RESIDUE POSITION      - KEEP FOR NOW, AS MAY RELATE TO INPUT TO MAKE EMBEDDINGS.
        'A_id',                     # ATOM POSITION         - SORT ON THIS, KEEP IN DF.
        'S_mon_id',                 # RESIDUE (3-letter)    - USE AS SANITY-CHECK AGAINST 'A_label_comp_id', KEEP IN DF.
        'A_label_comp_id',          # RESIDUE (3-letter)    - USE AS SANITY-CHECK WITH 'S_mon_id', THEN REMOVE IT.
        'A_label_atom_id',          # ATOM                  - KEEP IN DF.
        'A_label_asym_id',          # CHAIN                 - JUST REMOVE
        'S_asym_id',                # CHAIN                 - JUST KEEP
        'A_Cartn_x',                # COORDINATES           - X-COORDINATES
        'A_Cartn_y',                # COORDINATES           - Y-COORDINATES
        'A_Cartn_z',                # COORDINATES           - Z-COORDINATES
        'A_occupancy',              # OCCUPANCY             - FILTER ON THIS, THEN REMOVE IT.
        #'b_iso_or_equiv',           # B-FACTORS             - BUT THIS IS ONLY FOR CRYSTALLOGRAPHIC DATA NOT NMR
    ]]


def _join_atomsite_to_polyseq(atomsite: pd.DataFrame, polyseq: pd.DataFrame) -> pd.DataFrame:
    """
    pdf_merged MUST BE FOR ONE POLYPEPTIDE CHAIN ONLY.
    """
    return pd.merge(left=polyseq, right=atomsite, left_on=['S_seq_id'], right_on=['A_label_seq_id'], how='outer')


def _split_up_by_chain(atomsite_pdf: pd.DataFrame, polyseq_pdf: pd.DataFrame) -> list:
    """
    :return: List of tuples containing each given pdf for each polypeptide chain,
    e.g. [(`atomsite_pdf_A`, `polyseq_pdf_A`), (`atomsite_pdf_B`, `polyseq_pdf_B`, etc)]
    """
    chains = atomsite_pdf['A_label_asym_id'].unique()
    num_of_chains_A = len(chains)
    num_of_chains_S = polyseq_pdf['S_asym_id'].nunique()
    if num_of_chains_A != num_of_chains_S:
        print(f"There are {num_of_chains_A} chains in `_atom_site`, but {num_of_chains_S} chains in "
              f"`_pdbx_poly_seq_scheme`. Probably a chain with non-natural amino acids (which would have been"
              f" filtered out by preceding function `_remove_hetatm_rows`.")
    grouped_atomsite_dfs = [group_df for _, group_df in atomsite_pdf.groupby('A_label_asym_id')]
    grouped_polyseq_dfs = [group_df for _, group_df in polyseq_pdf.groupby('S_asym_id')]
    grouped_tuple = [(grp_as, grp_ps) for grp_as, grp_ps in zip(grouped_atomsite_dfs, grouped_polyseq_dfs)]
    assert len(chains) == len(grouped_tuple)
    return grouped_tuple


def _remove_hetatm_rows(atomsite_pdf: pd.DataFrame) -> pd.DataFrame:
    atomsite_pdf = atomsite_pdf.drop(atomsite_pdf[atomsite_pdf['A_group_PDB'] == 'HETATM'].index)
    # OR KEEP ONLY ROWS WITH 'ATOM' GROUP. NOT SURE IF ONE APPROACH IS BETTER THAN THE OTHER:
    # atom_site_pdf = atom_site_pdf[atom_site_pdf.A_group_PDB == 'ATOM']
    if atomsite_pdf.empty:
        print('Having just removed all HETATM rows, the df is now empty. '
              'This suggests this PDB did not contain any natural amino acids (e.g. Could be all Norleucine)')
    return atomsite_pdf


def extract_fields_from_atom_site(mmcif: dict) -> pd.DataFrame:
    _atom_site = '_atom_site.'
    group_pdbs = mmcif[_atom_site + 'group_PDB']                        # GROUP ('ATOM' or 'HETATM')
    ids = mmcif[_atom_site + 'id']                                      # ATOM POSITIONS
    label_atom_ids = mmcif[_atom_site + 'label_atom_id']                # ATOMS
    label_comp_ids = mmcif[_atom_site + 'label_comp_id']                # RESIDUE (3-LETTER)
    label_asym_ids = mmcif[_atom_site + 'label_asym_id']                # CHAIN
    label_seq_ids = mmcif[_atom_site + 'label_seq_id']                  # RESIDUE POSITION
    x_coords = mmcif[_atom_site + 'Cartn_x']                            # CARTESIAN X COORDS
    y_coords = mmcif[_atom_site + 'Cartn_y']                            # CARTESIAN Y COORDS
    z_coords = mmcif[_atom_site + 'Cartn_z']                            # CARTESIAN Z COORDS
    occupancies = mmcif[_atom_site + 'occupancy']                       # OCCUPANCY
    model_nums = mmcif[_atom_site + 'pdbx_PDB_model_num']               # MODEL NUMBER
    # b_factors = mmcif[_atom_site + 'b_iso_or_equiv']                    # B-FACTORS 
        
    # 'A_' IS FOR `_atom_site`
    atom_site = pd.DataFrame(
        data={
            'A_group_PDB': group_pdbs,                          # 'ATOM' or 'HETATM'
            'A_id': ids,                                        # 1,2,3,4,5,6,7,8,9,10, etc
            'A_label_atom_id': label_atom_ids,                  # 'N', 'CA', 'C', 'O', etc
            'A_label_comp_id': label_comp_ids,                  # 'ASP', 'ASP', 'ASP', etc
            'A_label_asym_id': label_asym_ids,                  # 'A', 'A', 'A', 'A', etc
            'A_label_seq_id': label_seq_ids,                    # 1,1,1,1,1,1,2,2,2,2,2, etc
            'A_Cartn_x': x_coords,                              # COORDS
            'A_Cartn_y': y_coords,                              # COORDS
            'A_Cartn_z': z_coords,                              # COORDS
            'A_occupancy': occupancies,                         # BETWEEN 0 AND 1.0
            'A_pdbx_PDB_model_num': model_nums,                 # BETWEEN 1 AND ANY NUMBER (TYPICALLY LESS THAN 40)
      #      'b_iso_or_equiv': b_factors                         # FLOATS 0 TO 100 ???
        })

    return atom_site


def extract_fields_from_poly_seq(mmcif: dict) -> pd.DataFrame:
    _pdbx_poly_seq_scheme = '_pdbx_poly_seq_scheme.'
    seq_ids = mmcif[_pdbx_poly_seq_scheme + 'seq_id']                               # RESIDUE POSITION
    mon_ids = mmcif[_pdbx_poly_seq_scheme + 'mon_id']                               # RESIDUE (3-LETTER)
    pdb_seq_nums = mmcif[_pdbx_poly_seq_scheme + 'pdb_seq_num']                     # RESIDUE POSITION
    asym_ids = mmcif[_pdbx_poly_seq_scheme + 'asym_id']                             # CHAIN

    # 'S_' IS FOR `_pdbx_poly_seq_scheme`
    poly_seq = pd.DataFrame(
        data={
            'S_seq_id': seq_ids,                                # 1,1,1,1,1,1,2,2,2,2,2, etc
            'S_mon_id': mon_ids,                                # 'ASP', 'ASP', 'ASP', etc
            'S_pdb_seq_num': pdb_seq_nums,                      # 1,1,1,1,1,1,2,2,2,2,2, etc
            'S_asym_id': asym_ids                               # 'A', 'A', 'A', 'A', etc
        })
    return poly_seq


def parse_pdb_snapshot(t):
    """
    Filter out non-alpha-carbon atoms and keep only `serial`, `x`, `y`, `z`, `resSeq`, `resName`.
    The PDB column names equivalents in mmCIFs:
    `serial`        <==> `_atom_site.id`
    `name`          <==> `_atom_site.label_atom_id`
    `x`, `y`, `z`   <==> `_atom_site.Cartn_x`, `_atom_site.Cartn_y`, `_atom_site.Cartn_z`
    `resSeq`        <==> `_pdbx_poly_seq_scheme.seq_id`
    `resName`       <==> `_pdbx_poly_seq_scheme.mon_id`
    :return:
    """
    pdf_chain = t.topology.to_dataframe()[0]
    coords = t.xyz[0]  # shape (n_atoms, 3)
    pdf_chain['x'] = coords[:, 0]
    pdf_chain['y'] = coords[:, 1]
    pdf_chain['z'] = coords[:, 2]
    pdf_chain = pdf_chain.loc[pdf_chain['name'].isin(('CA',))]
    pdf_chain = pdf_chain[['resSeq', 'resName', 'serial', 'x', 'y', 'z']]
    return pdf_chain


def parse_pdb_alpha_carbs_only(pdbid_chain: str, u):
    rpath_CA = os.path.join('..', 'data', 'ATLAS_parsed', pdbid_chain, 'CA_only')
    os.makedirs(rpath_CA, exist_ok=True)

    ca = u.select_atoms('protein and name CA')
    # Write all frames to separate files
    for i, _ in enumerate(u.trajectory):
        ca.write(os.path.join(rpath_CA, f'frame_{i:05d}.pdb'))


def parse_cif(pdb_id: str, mmcif_dict: dict) -> Tuple[List[pd.DataFrame], list]:
    print(f'parsing {pdb_id}')
    polyseq_pdf = extract_fields_from_poly_seq(mmcif_dict)
    atomsite_pdf = extract_fields_from_atom_site(mmcif_dict)
    atomsite_pdf = _remove_hetatm_rows(atomsite_pdf)
    # GENERATE A LIST OF TUPLES, EACH TUPLE IS THE ATOMSITE AND POLYSEQ DATA FOR A SINGLE CHAIN
    all_chains_pdfs = _split_up_by_chain(atomsite_pdf, polyseq_pdf)
    # IF CHAIN SPECIFIED IN PDBID, REMOVE ANY OTHER CHAINS <-- Implemented later in tokeniser.py, but could move here.
    parsed_cif_by_chain = []
    empty_pdbidchains = []
    for chain_pdf in all_chains_pdfs:
        atomsite_pdf, polyseq_pdf = chain_pdf
        joined_pdf_chain = _join_atomsite_to_polyseq(atomsite_pdf, polyseq_pdf)
        try:
            chain = joined_pdf_chain['S_asym_id'].iloc[0]
        except:
            print(f"joined_pdf_chain['S_asym_id'].iloc[0] on line 246 is failing for {pdb_id}")
            print('Leaving this one out.')
            print(f"joined_pdf_chain={joined_pdf_chain}")
            continue

        joined_pdf_chain = joined_pdf_chain.loc[joined_pdf_chain['A_label_atom_id'].isin(('CA',))]  # ALPHA-CARBON ONLY
        # print(f'joined_pdf_chain.shape after removing everything except alpha-carbons={joined_pdf_chain.shape}')
        if joined_pdf_chain.empty:
            pdbid_chain = f'{pdb_id}_{chain}'
            empty_pdbidchains.append(pdbid_chain)
            continue
        joined_pdf_chain = _rearrange_cols(joined_pdf_chain)
        joined_pdf_chain = _cast_number_strings_to_numeric_types(joined_pdf_chain)
        joined_pdf_chain = _cast_objects_to_stringdtype(joined_pdf_chain)
        joined_pdf_chain = _sort_by_chain_residues_atoms(joined_pdf_chain)
        joined_pdf_chain = _replace_low_occupancy_coords_with_nans(joined_pdf_chain)
        joined_pdf_chain = _process_missing_data(joined_pdf_chain, impute=False)

        joined_pdf_chain = joined_pdf_chain[['A_pdbx_PDB_model_num',                    # MODEL NUMBER
                                             'S_asym_id',                               # CHAIN
                                             'S_seq_id',                                # RESIDUE POSITION
                                             'S_mon_id',                                # RESIDUE NAME (3-LETTER)
                                             'A_id',                                    # ATOM POSITION
                                             'A_label_atom_id',                         # ATOM NAME
                                             'A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']]    # COORDINATES
                                            # 'b_iso_or_equiv']]                       # B-FACTORS (only for crystallographic data, not NMR).
        parsed_cif_by_chain.append(joined_pdf_chain)
    return parsed_cif_by_chain, empty_pdbidchains


# if __name__ == '__main__':
#     from Bio.PDB.MMCIF2Dict import MMCIF2Dict
#     pdbid = '2N2K'
#     cif_pdfs_per_chain = parse_cif(pdb_id=pdbid,
#                                    mmcif_dict=MMCIF2Dict(f'../../data/NMR/raw_cifs/heteromeric/{pdbid}.cif'))
#     pass