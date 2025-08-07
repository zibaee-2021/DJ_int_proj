import os, glob
import sys
from collections import defaultdict

def _rp_nmr_dir():
    return os.path.join('..', '..', 'data', 'NMR')



def _reduce_to_minimum_header(preamble: list) -> list:
    min_preamble = []
    lines2include = ['REMARK   2 RESOLUTION. NOT APPLICABLEx.']
    keepers = ['HEADER', 'TITLE', 'EXPDTA', 'AUTHOR', 'CRYST1']
    for line in preamble:
        linestart = line[:6].strip()
        linemodel = line[:15].strip()
        if (line not in lines2include and
                linestart not in keepers and
                linemodel != 'MODEL        1'):
            continue
        else:
            min_preamble.append(line)
    return min_preamble


def _further_parsing(rp_pdbchain_f: str):
    """
    write_pdb_file_per_chain() does not include MODEL .. # after ENDMDL for some chains (because models can include
    multiple chains that are not each separated out by ENDMDL\nMODEL.. # etc, but are required in my PDBchain files).
    So I am reading them each in, to add the 'MODEL ...#' as appropriate.
    """
    pdbc_lines = []
    previous_record = ''
    model_num = 1
    with open(rp_pdbchain_f, 'r') as f:
        lines = f.readlines()
        # is there a pyhthon function that will return the index podition of a given string, e.g.  'MODEL        1'
    with open(rp_pdbchain_f, 'r') as f:
        for i, line in enumerate(f):
            record = line[: 6].strip()
            if previous_record == 'ENDMDL' and record != 'MASTER':
                model_num += 1
                pdbc_lines.append(f'MODEL        {model_num}\n')
            else:
                pdbc_lines.append(line)
            previous_record = record

    with open(rp_pdbchain_f, 'w') as f:
        f.writelines(f'{line}' for line in pdbc_lines)
    print(f'Wrote: {os.path.basename(rp_pdbchain_f)}')


def _pdbid_dict_chain(pdbid_chains: list) -> dict:
    """
    Takes relative path of txt file that has list of sol NMR PDBid_chains with > 1 model.
    Makes a dict with unique PDBid as the key, mapped to all its chains, according to the list in the txt file.
    Returns a list of the PDBid-chains, and a dictionary mapping PDB ids to chains.
    """
    pdbid_chains_dict = defaultdict(list)
    pdbid_chains.sort()
    for pdbid_chain in pdbid_chains:
        pid, ch = pdbid_chain.split('_')
        pdbid_chains_dict[pid].append(ch)
    return dict(pdbid_chains_dict)



def write_pdb_file_per_chain(rp_pdb_f: str, chains_to_include: list, rp_dst_dir: str):
    """
    Splits a .pdb file into per-chain .pdb files.
    Each output file will contain all ATOM lines for a single chain,
    and will retain END/TER lines per model.
    """
    chain_lines = defaultdict(list)
    chain = ''
    preamble = list()
    previous_record = ''
    got_chain = False
    include_this_chain = True

    with open(rp_pdb_f, 'r') as f:
        all_lines = f.readlines()  # used for making preamble (~ line 101 below)

    with open(rp_pdb_f, 'r') as f:
        for i, line in enumerate(f):
            record = line[: 6].strip()
            if record in ['ATOM', 'TER', 'MODEL', 'ENDMDL']:
                if record == 'ATOM':
                    chain = line[21]
                    include_this_chain = chain in chains_to_include
                    got_chain = True
                if not include_this_chain:
                    continue
                if got_chain and previous_record.startswith('MODEL'):
                    if previous_record != 'MODEL        1':
                        chain_lines[chain].append(f'{previous_record}\n')
                if got_chain and record.startswith('ATOM') and previous_record.startswith('TER'):  # 1 model can >1 chain
                    if chain != previous_chain:
                        chain_lines[previous_chain].append(f'ENDMDL\n')
                if got_chain:
                    chain_lines[chain].append(line)
                if line[: 20].strip() == 'MODEL        1':
                    preamble = all_lines[0: i + 1]
                    preamble = _reduce_to_minimum_header(preamble)
            previous_record = line[: 20].strip()
            previous_chain = chain
    for chain in chain_lines.keys():
        chain_lines[chain] = preamble + chain_lines[chain]
        chain_lines[chain].append('MASTER\nEND\n')

    pdb_id = os.path.basename(rp_pdb_f).split('.')[0]
    for chain, lines in chain_lines.items():
        out_path = os.path.join(rp_dst_dir, f'{pdb_id}_{chain}.pdb')
        with open(out_path, 'w') as out_f:
            out_f.writelines(lines)
        print(f'Wrote: {os.path.basename(out_path)}')

if __name__ == '__main__':

    # with open(os.path.join(_rp_nmr_dir(), 'multimodel_lists', 'multimod_2713_hetallchains_hom1chain.lst'), 'r') as f:
    #     pidchains = f.read().splitlines()
    # pdbids = list(set(pidchain[:-2] for pidchain in pidchains))
    # pdbid_chains_dict = _pdbid_dict_chain(pidchains)
    #
    # rp_raw_pdbs_dir = os.path.join(_rp_nmr_dir(), 'raw_pdbs', 'hethom_combined')
    # rp_pdbchains_dst_dir = os.path.join(_rp_nmr_dir(), 'pdb_chains', 'multimod_2713_hetallchains_hom1chain')
    # os.makedirs(rp_pdbchains_dst_dir, exist_ok=True)
    #
    # for pdbid_, chains_to_include_ in pdbid_chains_dict.items():
    #     rp_raw_pidchain_f = os.path.join(rp_raw_pdbs_dir, f'{pdbid_}.pdb')
    #     write_pdb_file_per_chain(rp_pdb_f=rp_raw_pidchain_f, chains_to_include=chains_to_include_,
    #                              rp_dst_dir=rp_pdbchains_dst_dir)
    #     pass

    rp_pdbchains_dir = os.path.join(_rp_nmr_dir(), 'pdb_chains', 'multimod_2713_hetallchains_hom1chain')
    rp_pdb_files = sorted(glob.glob(os.path.join(rp_pdbchains_dir, '*.pdb')))
    for rp_pdb_f_ in rp_pdb_files:
        _further_parsing(rp_pdb_f_)
