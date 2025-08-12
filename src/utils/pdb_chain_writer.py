import os, glob
from typing import override
from pathlib import Path
from collections import defaultdict
import numpy as np
import MDAnalysis as mda
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select



def _rp_nmr_dir():
    return os.path.join('..', '..', 'data', 'NMR')

##### NOT USED NOW -- REPLACED BY BIOPYTHON -- ############################################
def __reduce_to_minimum_header(preamble: list) -> list:
    """
    Not used now -- REPLACED BY BIOPYTHON --
    """
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


def __further_parsing(rp_pdbchain_f: str):
    """
    Not used now -- REPLACED BY BIOPYTHON --
    write_pdb_file_per_chain() does not include MODEL .. # after ENDMDL for some chains (because models can include
    multiple chains that are not each separated out by ENDMDL\nMODEL.. # etc, but are required in my PDBchain files).
    So I am reading them each in, to add the 'MODEL ...#' as appropriate.
    """
    pdbc_lines = []
    previous_record = ''
    model_num = 1
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
    print(f'{os.path.basename(rp_pdbchain_f)} written.')


def __pdbid_dict_chain(pdbid_chains: list) -> dict:
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


def __write_pdb_file_per_chain(rp_pdb_f: str, chains_to_include: list, rp_dst_dir: str):
    """
    NOT USED NOW -- REPLACED BY BIOPYTHON --
    Splits a .pdb file into per-chain .pdb files.
    Each output file will contain all ATOM lines for a single chain,
    and will retain END/TER lines per model.
    """
    chain_lines = defaultdict(list)
    chain = ''
    preamble = list()
    previous_record = ''
    got_chain = False

    with open(rp_pdb_f, 'r') as f:
        all_lines = f.readlines()  # used for making preamble (~ line 101 below)

    with open(rp_pdb_f, 'r') as f:
        for i, line in enumerate(f):
            record = line[: 6].strip()
            if record in ['ATOM', 'TER', 'MODEL', 'ENDMDL']:
                if record == 'ATOM':
                    chain = line[21]
                    got_chain = True
                    pdb_id = os.path.basename(rp_pdb_f).removesuffix('.pdb')

                if got_chain and chain not in chains_to_include:
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
                    preamble = __reduce_to_minimum_header(preamble)
            previous_record = line[: 20].strip()
            previous_chain = chain
    for chain in chain_lines.keys():
        chain_lines[chain] = preamble + chain_lines[chain]
        chain_lines[chain].append('MASTER\nEND\n')

    pdb_id = os.path.basename(rp_pdb_f).split('.')[0]
    for chain, lines in chain_lines.items():

        if chain in chains_to_include:
            out_path = os.path.join(rp_dst_dir, f'{pdb_id}_{chain}.pdb')
            with open(out_path, 'w') as out_f:
                out_f.writelines(lines)
            print(f'Wrote: {os.path.basename(out_path)}')


def __write_pdb_hetatmout_caonly(rp_pdf_f: str, rp_dst_parsed_pdbs_dir: str) -> None:
    """
    NOT USED NOW -- REPLACED BY BIOPYTHON --
    """
    pidc = os.path.basename(rp_pdf_f).removeprefix('.pdb')
    print(pidc)
    u = mda.Universe(rp_pdf_f)  # can also be .cif
    ca_atoms = u.select_atoms('protein and name CA') # standard residues (no HETATM) and alpha-carbons only.

    rp_parsed_pdb_f = os.path.join(rp_dst_parsed_pdbs_dir, f'{pidc}')
    with mda.Writer(rp_parsed_pdb_f, ca_atoms.n_atoms) as writer:
        writer.write(ca_atoms)


def __write_mean_coords_to_pdb(rp_pdb_f: str) -> None:
    """
    NOT USED NOW -- REPLACED BY BIOPYTHON --
    Using MDAnalysis to calculate mean coords of given PDB/PDBchain file and write out to PDB/PDBchain file.
    Note: the PDB it writes out:
     - starts with 'MODEL 1' even though its an average of all the models.
     - lacks 'MASTER' at end and 'TER' at terminus.
    I assume it is still a valid PDB as it seems unlikely to me that this is an error by MDAnalysis.
    """
    u = mda.Universe(rp_pdb_f)
    n_atoms = len(u.atoms)
    coords_per_model = []  # store coords per model
    for _ in u.trajectory:  # iterate through each model
        coords_per_model.append(u.atoms.positions.copy())
    coords_array = np.array(coords_per_model)  # shape=(n_models, n_atoms, 3)
    mean_coords = coords_array.mean(axis=0)  # mean coords over all models
    u.atoms.positions = mean_coords  # set mean coords into "Universe"
    with mda.Writer(os.path.basename(rp_pdb_f), n_atoms) as W:
        W.write(u.atoms)


########################################################################################

def write_mean_models_to_pdb(rp_pidchain_f: str, dst_mean_coords_pdb_dir: str | None = None) -> None:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('avg', rp_pidchain_f)

    # Collect coords per (chain, full residue id, atom name)
    coord_map = {}  # key -> list of (3,) arrays
    for model in struct:  # iterate models
        for chain in model:
            for residue in chain:
                for atom in residue:
                    key = (chain.id, residue.id, atom.id)
                    coord_map.setdefault(key, []).append(atom.coord.copy())

    # Mean per key (average only over models where the atom exists)
    mean_coord = {k: np.vstack(v).mean(axis=0) for k, v in coord_map.items()}

    # Select that (a) keeps only model 0, (b) overwrites coords with the mean
    class MeanModel0(Select):
        def accept_model(self, model):
            return model.id == 0
        def accept_atom(self, atom):
            key = (atom.get_parent().get_parent().id,  # chain.id
                   atom.get_parent().id,               # residue.id tuple
                   atom.id)
            if key in mean_coord:
                atom.set_coord(mean_coord[key])
                return True
            return False  # atom missing from all models? skip it

    # Write out single-model, mean-coord PDB
    io = PDBIO()
    io.set_structure(struct)
    pidchain_pdbext = os.path.basename(rp_pidchain_f)
    rp_meancoords_pidchain_f = os.path.join(dst_mean_coords_pdb_dir, pidchain_pdbext)
    io.save(rp_meancoords_pidchain_f, select=MeanModel0())
    print(f'Wrote {pidchain_pdbext}')
    return dst_mean_coords_pdb_dir


class ChainAlphaCOnly(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    @override
    def accept_chain(self, chain):
        return chain.id == self.chain_id

    @override
    def accept_residue(self, residue):
        return residue.id[0] == ' '

    @override
    def accept_atom(self, atom):
        return atom.id == 'CA'


def parse_ca_atom_and_write_pdb_per_chain() -> None:
    # GET REL PATHS TO ALL 1724 RAW CIF FILES AVAILABLE IN THE TWO DATA DIRS (HET AND HOM):
    rp_rawcifs_het_dir = os.path.join(_rp_nmr_dir(), 'raw_cifs', 'heteromeric')
    rp_rawcifs_hom_dir = os.path.join(_rp_nmr_dir(), 'raw_cifs', 'homomeric')
    rp_rawcifs_het_1038_cifs = glob.glob(os.path.join(rp_rawcifs_het_dir, '*.cif'))
    rp_rawcifs_hom_686_cifs = glob.glob(os.path.join(rp_rawcifs_hom_dir, '*.cif'))
    rp_rawcifs_all_1724_cifs = rp_rawcifs_het_1038_cifs + rp_rawcifs_hom_686_cifs

    # GET ALL 2713 PID-CHAINS FROM LIST FILE:
    with open(os.path.join(_rp_nmr_dir(), 'multimodel_lists', 'multimod_2713_hetallchains_hom1chain.lst'), 'r') as f:
        mm_2713_pidchains = f.read().splitlines()
    mm_2713_pidchains.sort()
    length = len(mm_2713_pidchains)
    # DOUBLE-CHECK THEY ARE UNIQUE:
    length2 = len(set(mm_2713_pidchains))
    assert length == length2, 'duplicates in multimod_2713_hetallchains_hom1chain.lst ??'

    # EXTRACT UNIQUE PIDS FROM THE 2713 PID-CHAINS JUST TO CHECK IT MATCHES TO THOSE FROM _pdbid_dict_chain() BELOW:
    mm_1541_pids = list(set(pidchain[:-2] for pidchain in mm_2713_pidchains))
    mm_1541_pids.sort()
    length = len(mm_1541_pids)
    # BUILD DICT OF PIDS MAPPED TO LIST OF WANTED CHAINS. E.G. `['1A3C_A','1A3C_B',...]` TO `['1A3C': ['A','B'],...]`
    pidchains_dict = __pdbid_dict_chain(mm_2713_pidchains)
    length2 = len(pidchains_dict)
    assert length == length2, '_pdbid_dict_chain() not generating same number of unqiue pids as expected... ??'

    rp_pdbchains_dst_dir = os.path.join(_rp_nmr_dir(), 'pdb_chains', 'multimod_2713_hetallchains_hom1chain')
    os.makedirs(rp_pdbchains_dst_dir, exist_ok=True)

    for pdbid_, chains_to_include_ in pidchains_dict.items():
        print(f'{pdbid_} and {chains_to_include_} only.')
        parser = MMCIFParser(QUIET=True, auth_chains=False)
        rp_raw_cif = [rp_cif for rp_cif in rp_rawcifs_all_1724_cifs if pdbid_ == os.path.basename(rp_cif).removesuffix('.cif')][0]

        structure = parser.get_structure('struct', rp_raw_cif)

        # Get chain IDs from first model. (Note: these are not auth_asym_id, because I'm using mmCIFs, not legacy PDBs).
        first_model = next(structure.get_models())
        chain_ids = [chain.id for chain in first_model.get_chains()]

        # Save each chain in own PDB, keeping all models
        io = PDBIO()
        io.set_structure(structure)
        for cid in chain_ids:
            if cid in chains_to_include_:
                pidchain_pdbext = f'{pdbid_}_{cid}.pdb'
                rp_dst_pidchains_pdb_f = os.path.join(rp_pdbchains_dst_dir, pidchain_pdbext)
                io.save(rp_dst_pidchains_pdb_f, select=ChainAlphaCOnly(cid))
                print(f'Wrote {pidchain_pdbext}')


class OnlyModel(Select):
    def __init__(self, model_id):
        self.model_id = model_id  # Biopython's internal model id (usually 0..N-1)

    def accept_model(self, model):
        # Keep only the requested model
        return model.id == self.model_id

def split_models(rp_pidchain_f: str, rp_per_model_dst_dir: str | None = None) -> None:
    """
    Split a multi-model PDB into separate single-model PDB files.
    Returns list of output file paths.
    """
    rp_pidchain_f = Path(rp_pidchain_f)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('x', str(rp_pidchain_f))

    models = list(structure.get_models())
    assert len(models) >= 1, 'No models detected by Biopython'
    io = PDBIO()
    io.set_structure(structure)

    for i, model in enumerate(models, start=1):
        # i is 1-based for human-friendly filenames; model.id is Biopython's internal id
        rp_permodel_pidchain_dst_dir = os.path.join(rp_pidchains_per_model_dst_dir, rp_pidchain_f.stem)
        os.makedirs(rp_permodel_pidchain_dst_dir, exist_ok=True)
        pidchain_model_pdbext = f'{rp_pidchain_f.stem}_{i:02d}.pdb'
        rp_per_model_pidchain_f = os.path.join(rp_permodel_pidchain_dst_dir, pidchain_model_pdbext)
        io.save(str(rp_per_model_pidchain_f), select=OnlyModel(model_id=model.id))
        print(f'Wrote {pidchain_model_pdbext}')

    return None

if __name__ == '__main__':

    # parse_CA_ATOM_and_write_pdb_per_chain()
    # write_mean_models_to_pdb(rp_pidchain_f, rp_mean_coords_dst_dir)

    rp_pidchains_dir = os.path.join(_rp_nmr_dir(), 'pdb_chains', 'multimod_2713_hetallchains_hom1chain')

    rp_pidchains_per_model_dst_dir = os.path.join(rp_pidchains_dir, 'per_model')
    os.makedirs(rp_pidchains_per_model_dst_dir, exist_ok=True)

    rp_pidchain_files = glob.glob(os.path.join(rp_pidchains_dir, '*.pdb'))
    rp_pidchain_files.sort()

    for rp_pidchain_f_ in rp_pidchain_files:
        split_models(rp_pidchain_f_, rp_pidchains_per_model_dst_dir)