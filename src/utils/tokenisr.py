#!~/miniconda3/bin/python
"""
TOKENISER.PY
    - READ IN MMCIF FILE(S).
    - CALL `cif_parser.py` TO  PARSE MMCIF FILE(S)
    - TOKENISE MMCIF
    - WRITE TO .SSV FLATFILE.
----------------------------------------------------------------------------------------------------------------------
The following 14 mmCIF fields are extracted from the raw mmCIF files, parsed and tokenised into a dataframe.
The 14 fields are:

_atom_site:
    group_PDB           # 'ATOM' or 'HETATM'    - FILTER ON THIS, THEN REMOVE IT.
    label_seq_id        # RESIDUE POSITION      - JOIN TO S_seq_id, THEN REMOVE IT.
    label_comp_id       # RESIDUE (3-LETTER)    - FOR SANITY-CHECK WITH S_mon_id, THEN REMOVE IT.
    id                  # ATOM POSITION         - SORT ON THIS, KEEP IN DATAFRAME.
    label_atom_id       # ATOM                  - KEEP IN DATAFRAME.
    label_asym_id       # CHAIN                 - JOIN TO S_asym_id, KEEP IN DATAFRAME.
    Cartn_x             # COORDINATES           - ATOM X-COORDINATES
    Cartn_y             # COORDINATES           - ATOM Y-COORDINATES
    Cartn_z             # COORDINATES           - ATOM Z-COORDINATES
    occupancy           # OCCUPANCY             - FILTER ON THIS, THEN REMOVE IT.

_pdbx_poly_seq_scheme:
    seq_id              # RESIDUE POSITION      - JOIN TO A_label_seq_id. SORT ON THIS, KEEP IN DATAFRAME.
    mon_id              # RESIDUE (3-LETTER)    - USE FOR SANITY-CHECK AGAINST A_label_comp_id, KEEP IN DATAFRAME.
    pdb_seq_num         # RESIDUE POSITION      - KEEP FOR NOW, AS MAY RELATE TO INPUT TO MAKE EMBEDDINGS.
    asym_id             # CHAIN                 - JOIN TO A_label_asym_id, SORT ON THIS, THEN REMOVE IT.

----------------------------------------------------------------------------------------------------------------------
The output of the current `parse_tokenise_cif_write_flatfile()` function is a 12-column dataframe.
'_atom_site' is abbreviated to 'A_' prefix.
'_pdbx_poly_seq_scheme' is abbreviated to 'S_' prefix.
These 13 columns are:

A_label_asym_id       # CHAIN                 - JOIN ON THIS, SORT ON THIS, KEEP IN DF.
S_seq_id              # RESIDUE POSITION      - SORT ON THIS, KEEP IN DATAFRAME.
A_id                  # ATOM POSITION         - SORT ON THIS, KEEP IN DF.
A_label_atom_id       # ATOM                  - KEEP IN DF.
A_Cartn_x             # COORDINATES           - ATOM X-COORDINATES.
A_Cartn_y             # COORDINATES           - ATOM Y-COORDINATES.
A_Cartn_z             # COORDINATES           - ATOM Z-COORDINATES.
aa_label_num          # ENUMERATED RESIDUES   - EQUIVALENT TO `ntcodes` IN ORIGINAL RNA CODE. KEEP IN DF.
bb_or_sc              # BACKBONE OR SIDE-CHAIN ATOM ('bb' or 'sc'), KEEP FOR POSSIBLE SUBSEQUENT OPERATIONS.
bb_atom_pos           # ATOM POSITION CA OR MOST C-TERM OTHER BB ATOM, PER RESIDUE. KEEP IN DF.
bbindices             # INDEX POSITION OF THE ATOM POSITION (`A_id`) OF ALLOCATED BACKBONE ATOM.
atom_label_num        # ENUMERATED ATOMS      - EQUIVALENT TO `atomcodes` IN ORIGINAL RNA CODE. KEEP IN DF.
"""


import os
import re
import glob
import json
from typing import List, Tuple
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def get_nums_of_missing_data(pdf):
    pd_isna = int((pdf.map(lambda x: isinstance(x, float) and pd.isna(x))).sum().sum())
    pd_na = int((pdf.map(lambda x: x is pd.NA)).sum().sum())
    pd_nat = int((pdf.map(lambda x: x is pd.NaT)).sum().sum())
    np_nan = int((pdf.map(lambda x: x is np.nan)).sum().sum())
    none = int((pdf.map(lambda x: x is None)).sum().sum())
    empty_str = int((pdf.map(lambda x: x == ' ')).sum().sum())
    na_str = int((pdf.map(lambda x: x == 'na' or x == 'NA')).sum().sum())
    nan_str = int((pdf.map(lambda x: x == 'nan' or x == 'NAN' or x == 'NaN')).sum().sum())
    none_str = int((pdf.map(lambda x: x == 'none' or x == 'None')).sum().sum())
    null_str = int((pdf.map(lambda x: x == 'null' or x == 'Null')).sum().sum())

    counts = {
        'NaN': pd_isna,
        'pd.NA': pd_na,
        'pd.NaT': pd_nat,
        'np.nan': np_nan,
        'None': none,
        ' ': empty_str,
        "na": na_str,
        "nan": nan_str,
        "none": none_str,
        "null": null_str
    }
    print(counts)
    has_missing_data = any(value > 0 for value in counts.values())

    if has_missing_data:
        nan_positions = pdf.map(lambda x: (isinstance(x, float) and pd.isna(x) or
                                           x is pd.NA or
                                           x is pd.NaT or
                                           x is np.nan or
                                           x is None or
                                           x == ' ' or
                                           x in ['na', 'NA'] or
                                           x in ['nan', 'NAN', 'NaN'] or
                                           x in ['none', 'None'] or
                                           x in ['null', 'Null']
                                           ))

        # Identify the row (index) and column names where the above checks are True
        nan_positions = nan_positions.stack().reset_index()  # Convert to long format
        nan_positions = nan_positions[nan_positions[0]]  # Filter only rows where True
        nan_positions = nan_positions[['level_0', 'level_1']]  # Keep only index and column info
        nan_positions.columns = ['index', 'column']  # Rename columns
        raise ValueError('There are missing values.. needs to be addressed.')

    # missing_strings = ['NaN', 'None', 'N/A', 'missing', 'NULL', '']
    # pdf = pdf.replace(missing_strings, np.nan)
    return pdf


def _each_column_has_expected_values(pdf_chain):
    # THIS IS AN EXTRA VALIDATION STEP THAT IS NOT NEEDED FOR THE 565 PROTEINS, BUT WOULD BE USEFUL FOR FUTURE
    # WILL CHECK THAT EACH NUMERIC COLUMN HAS VALUES IN THE EXPECTED RANGE.
    pass


def _chdir_to_data_layer() -> str:
    """
    Change current working dir to `data_layer`. (Intended to be switched back at end of calling function).
    :return: Current working directory *before* it had been changed to local dir (`data_handler`), to enable calling
    function to restore it.
    """
    cwd = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    return cwd


def _read_json_from_data_dir(fname: str) -> dict:
    """
    Read given json file from diffSock/data/{fname} to a Python dict.
    :param fname: File name of json file to read. IMPORTANT: Subdir paths are expected to be included.
    e.g. 'enumeration/fname.json' without starting forward slash.
    :return: The read-in json file, as a Python dict.
    """
    cwd = _chdir_to_data_layer()
    fname = fname.removeprefix('/').removesuffix('.json')
    relpath_json = f'../data/{fname}.json'
    assert os.path.exists(relpath_json)
    try:
        with open(relpath_json, 'r') as json_f:
            my_dict = json.load(json_f)
    except FileNotFoundError:
        print(f'{my_dict} does not exist.')
    except Exception as e:
        print(f"An error occurred: {e}")

    _restore_original_working_dir(cwd)
    return my_dict


def _read_enumerations_json(fname: str) -> dict:
    fname = fname.removesuffix('.json')
    return _read_json_from_data_dir(fname=f'enumeration/{fname}.json')


def _enumerate_atoms(pdf: pd.DataFrame) -> pd.DataFrame:
    """
    Enumerate atoms of mmCIF for one protein, by mapping via pre-written json at `data/enumeration`. Add
    this enumeration to a new column `atom_label_num`. It serves as the tokenised form of polypeptide atoms for this
    protein, to be read later to `atomcodes` array.
    :param pdf: Dataframe of one protein mmCIF, containing atoms to enumerate to new column.
    :return: Given dataframe with one new column holding the enumerated atoms data. Expected to have 12 columns.
    """
    # MAKE NEW COLUMN FOR ENUMERATED ATOMS ('C', 'CA', ETC), USING JSON->DICT, CAST TO INT:
    atoms_enumerated = _read_enumerations_json(fname='unique_atoms_only_no_hydrogens')
    pdf['atom_label_num'] = (pdf['A_label_atom_id']
                             .map(atoms_enumerated)
                             .astype('Int64'))
    expected_num_of_cols = 12
    assert len(pdf.columns) == expected_num_of_cols, \
        f'Dataframe should have {expected_num_of_cols} columns. But this has {len(pdf.columns)}'
    return pdf


def _enumerate_residues(pdf: pd.DataFrame) -> pd.DataFrame:
    # MAKE NEW COLUMN FOR ENUMERATED RESIDUES, USING JSON->DICT, CAST TO INT.
    # `residues_enumerated` DICT KEY AND `S_mon_id` COLUMN VALUES MAP VIA 3-LETTER RESIDUE NAMES:
    residues_enumerated = _read_enumerations_json(fname='residues')
    pdf.loc[:, 'aa_label_num'] = (pdf['S_mon_id'].map(residues_enumerated).astype('Int64'))
    expected_num_of_cols = 11
    assert len(pdf.columns) == expected_num_of_cols, \
        f'Dataframe should have {expected_num_of_cols} columns. But this has {len(pdf.columns)}'
    return pdf


def enumerate_atoms_and_residues(pdfs: List[pd.DataFrame], pdb_id: str) -> List[pd.DataFrame]:
    """
    Enumerate residues, atoms, and residue-atoms pairs of given protein mmCIF data, and store to new columns in given
    dataframe. Currently hard-coded to only use atom data that lacks all hydrogen atoms.
    :param pdfs: List of dataframes of one protein mmCIF, per chain, containing atoms to enumerate to new columns.
    :return: Dataframe with new columns of enumerated data for residues, atoms, and residue-atoms pairs.
    """
    if pdb_id:
        print(f'PDBid={pdb_id}: enumerate atoms and residues')
    result_pdfs = list()
    for pdf in pdfs:
        pdf = _enumerate_residues(pdf)
        pdf = _enumerate_atoms(pdf)
        # pdf = _enumerate_residues_atoms(pdf)  # CURRENTLY NOT BEING USED.
        result_pdfs.append(pdf)
    return result_pdfs


def select_chains_to_use(pdfs: List[pd.DataFrame], chains: list = None, pdb_id: str = None) -> List[pd.DataFrame]:
    """
    Select which chains to keep for further parsing and tokenisation and to be written to flatfile. If no chains
    specified, all protein chains will be kept. (If no PDBid is given, just don't print anything).
    NOTE: CURRENTLY, THIS FUNCTION IS ONLY PERFORMING A SIMPLE CHAIN SELECTION LOGIC, NAMELY TO KEEP THE FIRST IN THE
    GIVEN LIST. FUTURE IMPLEMENTATION MAY LOOK AT KEEPING CHAINS THAT HAVE LESS THAN SOME GIVEN (E.G. 30 %) IDENTITY.
    :param pdfs: Dataframes, one per chain, for given protein mmCIF.
    :param chains: One of more chain(s) to keep. e.g. [A, C]. <Currently not used>
    :param pdb_id: Only for printing out as part of tracking function calls.
    :return: A list of one or more dataframes, according to which chains to keep (This is also temporary to avoid
    breaking subsequent operations).
    """
    if pdb_id:
        print(f'PDBid={pdb_id}: select which chains to keep.')
    if len(pdfs) > 0:
        pdfs = [pdfs[0]]
    return pdfs


def only_keep_chains_with_enough_bb_atoms(pdfs: List[pd.DataFrame], pdb_id: str = None) -> List[pd.DataFrame]:
    """
    Remove chain(s) from list of chains for this mmCIF if does not have more than a specified number of backbone
    polypeptides.
    :param pdfs: List of dataframes, one per chain, for one mmCIF.
    :param pdb_id: The PDBid of the corresponding mmCIF data.
    :return: List of dataframes, one per chain of protein, with any chains removed if they have less than the minimum
    permitted ratio of missing atoms (arbitrarily chosen for now).
    """
    # If more than this proportion of residues have no backbone atoms, remove the chain.
    MIN_RATIO_MISSING_BACKBONE_ATOMS = 0.0

    if pdb_id:
        print(f'PDBid={pdb_id}: remove chains without enough backbone atoms')
    result_pdfs = []
    for pdf in pdfs:
        aa_w_no_bbatom_count = (pdf
                                .groupby('S_seq_id')['bb_or_sc']
                                .apply(lambda x: 'bb' not in x.values)
                                .sum())
        total_atom_count = pdf.shape[0]
        if (aa_w_no_bbatom_count / total_atom_count) <= MIN_RATIO_MISSING_BACKBONE_ATOMS:
            result_pdfs.append(pdf)
    if len(result_pdfs) == 0:
        print(f'PDBid={pdb_id}: After removing chains that have too many residues that lack any backbone atom '
              f'coordinates at all, there are no chains left for this protein, PDBid={pdb_id}. It needs to be removed '
              f'from the dataset entirely.')
    return result_pdfs


def only_keep_chains_of_polypeptide(pdfs: List[pd.DataFrame], pdb_id: str) -> List[pd.DataFrame]:
    """
    Remove chain(s) from list of chains for this mmCIF if not polypeptide. Takes advantage of output of previous function
    which assigns to a new column called `bb_or_sc` 'bb' for polypeptide backbone atoms, 'sc' for polypeptide
    side-chain atoms, otherwise `pd.NaT` which therefore indicates the chain is not polypeptide.
    (RCSB mmCIF assigns a different chain to different molecule in a complex, e.g. RNA-protein complex, so if one row in
    the aforementioned column is `pd.NaT`, you should find that all are `pd.NaT`).
    :param pdfs: List of dataframes, one per chain, for one mmCIF.
    :param pdb_id: The PDBid of the corresponding mmCIF data.
    :return: List of dataframes, one per chain of protein, with any non-polypeptide chains removed.
    """
    if pdb_id:
        print(f'PDBid={pdb_id}: remove non-protein chains')
    result_pdfs = []
    for pdf in pdfs:
        try:
            chain = pdf['S_asym_id'].iloc[0]
            atleast_one_row_isna = pdf['bb_or_sc'].isna().any()
            all_rows_isna = pdf['bb_or_sc'].isna().all()
            if atleast_one_row_isna:
                nat_indices = pdf[pd.isna(pdf['bb_or_sc'])].index
                print(f'nat_indices={nat_indices}')
                print(f'There are atoms in chain={chain} of PDBid={pdb_id} which are not polypeptide atoms, so '
                      f'this chain will be excluded.')
                if not all_rows_isna:
                    print(f"It seems that while at least one row in column 'bb_or_sc' has null, "
                          f'not all rows are null. This is unexpected and should be investigated further. '
                          f'(Chain {chain} of PDBid {pdb_id}).')
            else:
                result_pdfs.append(pdf)
            if len(result_pdfs) == 0:
                print(f'PDBid={pdb_id}: After removing all non-polypeptide chains, there are no chains left. '
                      f'This should not occur, so needs to be investigated further.')
        except IndexError:
            print(f" `chain = pdf['S_asym_id'].iloc[0]` fails. "
                  f"\npdf['S_asym_id']={pdf['S_asym_id']}")
    return result_pdfs


def only_keep_alpha_carbons(pdfs: List[pd.DataFrame]) -> List[pd.DataFrame]:
    ALPHA_CARBON = ('CA',)
    result_pdfs = []

    for pdf in pdfs:
        _pdf = pdf.loc[pdf['A_label_atom_id'].isin(ALPHA_CARBON)]
        result_pdfs.append(_pdf)
    return result_pdfs


def make_new_bckbone_or_sdchain_col(pdfs: List[pd.DataFrame], pdb_id: str=None) -> List[pd.DataFrame]:
    BACKBONE = ('N', 'CA', 'C', 'O', 'OXT')  # AKA "MAIN-CHAIN"
    SIDECHAIN = ('CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CH2', 'CZ', 'CZ2',
                 'CZ3', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OD1', 'OD2', 'OE1', 'OE2', 'OG',
                 'OG1', 'OG2', 'OH', 'SD', 'SG')

    if pdb_id:
        print(f'PDBid={pdb_id}: make new column `bb_or_sc` - indicates whether atom is backbone or side-chain.')

    result_pdfs = list()

    for pdf in pdfs:
        is_backbone_atom = pdf['A_label_atom_id'].isin(BACKBONE)
        is_sidechain_atom = pdf['A_label_atom_id'].isin(SIDECHAIN)

        # MAKE NEW COLUMN TO INDICATE IF ATOM IS FROM POLYPEPTIDE BACKBONE ('bb) OR SIDE-CHAIN ('sc'):
        pdf.loc[:, 'bb_or_sc'] = np.select([is_backbone_atom, is_sidechain_atom], ['bb', 'sc'], default='placeholder')
        pdf.loc[pdf['bb_or_sc'] == 'placeholder', 'bb_or_sc'] = pd.NA

        expected_num_of_cols = 9
        assert len(pdf.columns) == expected_num_of_cols, \
            f'Dataframe should have {expected_num_of_cols} columns. But this has {len(pdf.columns)}'
        result_pdfs.append(pdf)

    return result_pdfs


def __split_pdbid_chain(pdbid_chain):
    match = re.match(r"^(.*)_([A-Za-z])$", pdbid_chain)
    if match:
        pdbid, chain = match.groups()
        return pdbid, chain
    else:
        return pdbid_chain, None


def _is_already_tokenised(relpath_toknsd_dir: str, pdbid: str) -> Tuple[bool, str]:
    """
    Check if pre-tokenised .ssv files of given PDBid/PDBid_chain are in given tokenised dir.
    :param relpath_toknsd_dir: Relative path to tokenised dir to look in.
    :param pdbid: PDBid or PDBid_chain to search for.
    :return: True if already tokenised .ssv file found in specified tokenised dir. Return PDBid_chain if True.
    """
    # MAKE A LIST OF PDBIDS THAT ARE SSVS IN TOKENISED DIR (I.E. HAVE ALREADY BEEN TOKENISED):
    pdbid_chain_pretokenised_list = []
    if relpath_toknsd_dir is not None:
        if os.path.exists(relpath_toknsd_dir):
            pdbid_chain_pretokenised_list = [item.removesuffix('.ssv')
                                            for item in os.listdir(relpath_toknsd_dir)
                                            if item.endswith('.ssv')]

    pdbid_chain_ptkn = ''
    is_pretokenised = False

    if pdbid in pdbid_chain_pretokenised_list:
        is_pretokenised = True
        pdbid_chain_ptkn = pdbid
    else:
        # IF PRE-TOKENISED AND SAVED AS ssv, PDBID WILL INCLUDE CHAIN SUFFIX.
        # But PDBids passed in may lack the chain, so make dict of PDBid to PDBid_chain for tokenised ones only.
        pdbid_to_pdbidchain_prtknsd = {}
        for pdbid_chain_prtknsd in pdbid_chain_pretokenised_list:
            pdbid_only, _ = __split_pdbid_chain(pdbid_chain_prtknsd)
            pdbid_to_pdbidchain_prtknsd[pdbid_only] = pdbid_chain_prtknsd
        pdbid_only, pdbid_chain_prtknsd = None, None
        # if PDBid only wo chain, check if in pretokenised list, after pdbid_chain str has chain part removed.
        if pdbid in pdbid_to_pdbidchain_prtknsd:
            is_pretokenised = True
            pdbid_chain_ptkn = pdbid_to_pdbidchain_prtknsd[pdbid]

    if is_pretokenised:
        print(f'PDBid = {pdbid} is already tokenised.')
    return is_pretokenised, pdbid_chain_ptkn


def _generate_list_of_pdbids_in_cif_dir(path_cif_dir: str) -> list:
    cifs = glob.glob(os.path.join(path_cif_dir, f'*.cif'))
    path_cifs = [cif.upper() for cif in cifs if os.path.isfile(cif)]
    pdb_id_list = []

    for path_cif in path_cifs:
        cif_basename = os.path.basename(path_cif)
        pdbid = os.path.splitext(cif_basename)[0]
        pdb_id_list.append(pdbid)
    return pdb_id_list


def _chdir_to_preprocessing_funcs_dir() -> None:
    """
    Change current working dir to `preprocessing_funcs`. .
    """
    cwd = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print(f'Changed working dir to {os.getcwd()} from {cwd}')


def _restore_original_working_dir(wd: str) -> None:
    """
    Change working directory to get working directory. Intended to work in tandem with `_chdir_to_data_layer()`.
    :param wd: Working directory to change to.
    """
    os.chdir(wd)


def write_tokenised_cifs_to_flatfiles(pdb_id: str, pdfs: List[pd.DataFrame], dst_data_dir=None, flatfile_format=None):
    """
    Write dataframe of single protein, and single chain, with columns CIF.S_seq_id, CIF.S_mon_id, CIF.A_id,
    CIF.A_label_atom_id, ColNames.MEAN_CORR_X, ColNames.MEAN_CORR_Y, ColNames.MEAN_CORR_Z to flat file(s) in local
    relative dir 'src/diffusion/diff_data/tokenised/' or to the top-level `data` dir.
    :param pdb_id: PDB id.
    :param pdfs: List of dataframes to write to flat file(s). One dataframe per polypeptide chain. TODO: hacked to one.
    :param dst_data_dir: Relative path to destination dir of flatfile of tokenised cif. (Will be called from either
    `diffSock/test`, `diffSock/src/preprocessing_funcs` or `diffSock/src/diffusion`. The responsibility for determining
    the relative destination path is left to the caller).
    :param flatfile_format: Flat file format to write to. E.g. 'ssv', 'csv' or 'tsv'.
    """
    flatfile_format = flatfile_format.removeprefix('.').lower()
    print(f'PDBid={pdb_id}: write tokenised to flatfile')
    for pdf in pdfs:
        if flatfile_format is None:
            flatfile_format = ['ssv']
        elif isinstance(flatfile_format, str):
            flatfile_format = [flatfile_format]
        cwd = ''  # to return to at end of this function.
        if not dst_data_dir:  # i.e. use the top-level general-use `data` dir & define relpath from data_layer
            print(f'You did not pass any destination dir path for writing the tokenised cif flat flatfile to. '
                  f'Therefore it will be written to the top-level data dir (`diffSock/data/tokenised`).')
            cwd = _chdir_to_data_layer()
            dst_data_dir = '../data/tokenised'
        else:
            os.makedirs(dst_data_dir, exist_ok=True)

        chain = pdf['S_asym_id'].unique()
        chain = chain[0]
        for flatfile in flatfile_format:
            sep = ' '
            if flatfile == 'tsv':
                sep = '\t'
            elif flatfile == 'csv':
                sep = ','
            dst_data_dir = dst_data_dir.removesuffix('/')
            pdb_id = pdb_id.removesuffix('.cif')
            pdf.to_csv(path_or_buf=f'{dst_data_dir}/{pdb_id}_{chain}.{flatfile}', sep=sep, index=False)

        # # For a more human-readable set of column-names:
        # pdf_easy_read = pdf.rename(columns={'S_seq_id': 'SEQ_ID',
        #                                     'S_mon_id': 'RESIDUES',
        #                                     'A_id': 'ATOM_ID',
        #                                     'A_label_atom_id': 'ATOMS',
        #                                     'MEAN_CORR_X': 'X',
        #                                     'MEAN_CORR_Y': 'Y',
        #                                     'MEAN_CORR_Z': 'Z'})
        # pdf_easy_read.to_csv(path_or_buf=f'../data/tokenised/easyRead_{pdb_id}.tsv', sep='\t',
        # index=False, na_rep='null')

        if not dst_data_dir:
            _restore_original_working_dir(cwd)


