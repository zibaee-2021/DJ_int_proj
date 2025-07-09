import glob
from time import time
import pandas as pd
from Bio.PDB import MMCIFParser, MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
parser = MMCIFParser(QUIET=True)


def generate_stats(cif_files):
    start = time()
    parser = MMCIFParser(QUIET=True)
    data = []
    for cif_file in cif_files:
        # Get number of models from the parsed structure:
        structure = parser.get_structure('', cif_file)
        num_models = len(list(structure))
        # Parse for metadata (publication year):
        cif_dict = MMCIF2Dict.MMCIF2Dict(cif_file)
        year = cif_dict.get('_citation.year', ['NA'])[0]

        chains_set = set()
        protein_chains_set = set()

        # Total number of all chains and total number of protein chains (using first model only):
        first_model = next(structure.get_models())
        for chain in first_model:
            chains_set.add(chain.id)
            is_protein = False
            for residue in chain:
            #     if residue.id[0] == ' ':  # standard residue
                if is_aa(residue, standard=True):
                    print(f"Chain {chain.id}, Residue {residue.get_resname()} is a standard amino acid")
                    is_protein = True
                    break
            if is_protein:
                protein_chains_set.add(chain.id)

        # Count number of alpha-carbons in each chain:
        for chain in first_model:
            ca_count = sum(1 for atom in chain.get_atoms() if atom.get_name() == 'CA')
            data.append({
                'File': cif_file,
                'Chain ID': chain.id,
                '# Alpha Carbons': ca_count,
                '# Chains Total': len(chains_set),
                '# Protein Chains': len(protein_chains_set),
                '# Models': num_models,
                'Year': year
            })
    df = pd.DataFrame(data)
    print(df)
    print(f'Completed {len(cif_files)} CIFs in {round(time() - start)} seconds')
    return df


if __name__ == '__main__':
    cifs = glob.glob('../data/NMR/raw_cifs/heteromeric/*.cif')
    # cifs = glob.glob('../data/NMR/raw_cifs/heteromeric/*.cif')
    cifs.sort()
    pdf = generate_stats(cif_files=cifs)
    pass