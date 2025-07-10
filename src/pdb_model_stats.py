import os, glob
from time import time
import pandas as pd
from Bio.PDB import MMCIFParser, MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
parser = MMCIFParser(QUIET=True)


def generate_stats(rpath_cifs):
    start = time()
    parser = MMCIFParser(QUIET=True)
    data = []
    for rpath_cif in rpath_cifs:
        structure = parser.get_structure('', rpath_cif)
        num_models = len(list(structure))

        cif_dict = MMCIF2Dict.MMCIF2Dict(rpath_cif)
        year = cif_dict.get('_citation.year', ['NA'])[0]

        chains_set = set()
        protein_chains_set = set()

        first_model = next(structure.get_models())
        for chain in first_model:
            chains_set.add(chain.id)
            is_protein = False
            for residue in chain:
                if is_aa(residue, standard=True):
                    # print(f'Chain {chain.id}, Residue {residue.get_resname()} is a standard amino acid')
                    is_protein = True
                    break
            if is_protein:
                protein_chains_set.add(chain.id)

        for chain in first_model:
            ca_count = sum(1 for atom in chain.get_atoms() if atom.get_name() == 'CA')
            data.append({
                'PDBid': os.path.basename(rpath_cif).removesuffix('.cif'),
                'Chain': chain.id,
                'CA count': ca_count,
                'Chain count': int(len(chains_set)),
                'Protein chain count': int(len(protein_chains_set)),
                'Model count': num_models,
                'Year': int(year)
            })
    pdf = pd.DataFrame(data)
    pdf_sorted = pdf.sort_values(by=['Year', 'Model count'], ascending=[True, True])
    print(pdf_sorted.dtypes)
    for col_to_cast in ['PDBid', 'Chain']:
        pdf_sorted[col_to_cast] = pdf_sorted[col_to_cast].astype("string")
    print(pdf_sorted.dtypes)
    print(f'Completed {len(rpath_cifs)} CIFs in {round(time() - start)} seconds')
    return pdf_sorted


if __name__ == '__main__':
    cifs = glob.glob('../data/NMR/raw_cifs/heteromeric/*.cif')
    # cifs = glob.glob('../data/NMR/raw_cifs/heteromeric/*.cif')
    cifs.sort()
    pdf = generate_stats(rpath_cifs=cifs[:4])
    print(pdf)
    pass