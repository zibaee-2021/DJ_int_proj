import os, glob
from time import time
import pandas as pd
from Bio.PDB import MMCIFParser, MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
parser = MMCIFParser(QUIET=True)
import mmseqs2

def generate_stats(rp_pdbid_chains: str, meric: str):
    start = time()
    parser = MMCIFParser(QUIET=True)
    data = []
    with open(rp_pdbid_chains, 'r') as f:
        pdbid_chains = [line.rstrip('\n') for line in f]
    pdbid_chains.sort()
    for pdbid_chain in pdbid_chains:
        pdbid, chain = pdbid_chain.split('_')
        rpath_cif = f'../data/NMR/raw_cifs/{meric}/{pdbid}.cif'
        structure = parser.get_structure('', rpath_cif)
        num_models = len(list(structure))
        year = structure.header['deposition_date'][:4]

        # cif_dict = MMCIF2Dict.MMCIF2Dict(rpath_cif)  # a few return '?', so I use 'get_structure()' despite some differences.
        # year2 = cif_dict.get('_citation.year', ['NA'])[0]
        # print(f'deposition={year}; citation={year2}')

        chains_set = set()
        protein_chains_set = set()

        first_model = next(structure.get_models())  # assumes first model can represent all models for this pdbid_chain
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

            if year != 'NA' and year != '?':
                print(f'year is given as {year}. Replacing with 0')
                year = int(year)
            else:
                year = 0

            chainpairs_mmseqs_ = mmseqs2.calc_mmseq(pdbid, chain)

            data.append({
                'PDBid': os.path.basename(rpath_cif).removesuffix('.cif'),
                'het_hom': os.path.basename(os.path.dirname(rpath_cif))[:3],
                'year': year,
                '#model': num_models,
                '#chain': int(len(chains_set)),
                '#prot_chain': int(len(protein_chains_set)),
                'chain': chain.id,
                'CA_count': ca_count,
                'identity': 89,  # TODO
                'similarity': 100  # TODO
            })
    pdf = pd.DataFrame(data)
    pdf_sorted = pdf.sort_values(by=['year', '#model'], ascending=[True, True])
    print(pdf_sorted.dtypes)
    for col_to_cast in ['PDBid', 'chain']:
        pdf_sorted[col_to_cast] = pdf_sorted[col_to_cast].astype('string')
    print(pdf_sorted.dtypes)
    print(f'Completed {len(pdbid_chains)} CIFs in {round(time() - start)} seconds')
    return pdf_sorted


if __name__ == '__main__':
    _meric = 'heteromeric'
    rpath_pdbid_chains = (f'../data/NMR/multimodel_PDBids/het_multimod_2104_pdbid_chains.txt')
    generate_stats(rp_pdbid_chains=rpath_pdbid_chains, meric=_meric)
    # pdf = generate_stats(rp_pdbid_chains=rpath_pdbid_chains, meric=_meric)

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

    pass