import os, glob
from time import time
from collections import defaultdict
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict
import mmseqs2

def generate_stats(rp_pdbid_chains: str=None, use_mmcif=True):
    """
    After parsing raw mmCIFs, the resulting PDBid_chains were written to a .lst file, which is used as follows:

    1. Read raw mmCIFs via Bio.PDB.MMCIFParser for each PDBid in .lst to get:
        - deposition year for each PDBid
        - total number of chains for each PDBid.

    2. Read PDBid_chains .lst file to get:
        - number of protein chains per PDBid.

    3. Read parsed SSV files to dataframes for each PDBid_chain in .lst to get:
        - number of models per PDBid_chain.
        - number of alpha carbons per PDBid_chain.

    4. Calculate sequence identity between every chain for each PDBid from list, by MMseqs2.
    """
    start = time()
    if not use_mmcif:
        parser = PDBParser(QUIET=True)
    else:
        parser = MMCIFParser(QUIET=True)
    stats, pdbids = [], []

    if not rp_pdbid_chains:
        rp_pdbid_chains = os.path.join('..', 'data', 'NMR', 'multimodel_PDBids', 'het_multimod_2104_pdbid_chains.txt')
    het_hom = os.path.basename(rp_pdbid_chains)[:3]
    het_hom = 'heteromeric' if het_hom == 'het' else 'homomeric'

    with open(rp_pdbid_chains, 'r') as f:
        pdbid_chains = [line.rstrip('\n') for line in f]
    pdbid_chains.sort()  # expected to already be sorted (i.e. before writing out), but to avoid relying on other code.
    pdbid_chains_dict = defaultdict(list)

    for pdbid_chain in pdbid_chains:
        pdbid, chain = pdbid_chain.split('_')
        pdbids.append(pdbid)
        pdbid_chains_dict[pdbid].append(chain)
    # 2104 heteromeric pdbid_chains --> 932 heteromeric pdbids.
    pdbids = list(set(pdbids))
    pdbids.sort()
    pdbid_chains_dict = dict(pdbid_chains_dict)

    rp_raw_cif_dir = os.path.join('..', 'data', 'NMR', 'raw_cifs', het_hom)
    rp_token_cif_dir = os.path.join('..', 'data', 'NMR', 'tokenised_cifs', het_hom)
    for pdbid in pdbids:
        rp_cif = os.path.join(rp_raw_cif_dir, f'{pdbid}.cif')
        bio_struct = parser.get_structure('', rp_cif)
        total_chain_count = len(list(next(bio_struct.get_models()).get_chains()))
        year = bio_struct.header['deposition_date'][:4]

        rp_fasta_file = mmseqs2.write_fasta(het_hom, pdbid, pdbid_chains_dict, bio_struct)
        result_pdf = mmseqs2.run_mmseqs_all_vs_all(rp_fasta_file, pdbid)
        filtered_pdf = mmseqs2.filter_results(result_pdf)
        identity = 0
        if not filtered_pdf.empty:
            pass
            # TODO generate dict of chain key to non-zero identity pairings as value.

        for chain in pdbid_chains_dict[pdbid]:
            pdbid_chain = f'{pdbid}_{chain}'
            rp_tok_cif = os.path.join(rp_token_cif_dir, f'{pdbid_chain}.ssv')
            pdbid_chain_pdf = pd.read_csv(rp_tok_cif, sep=' ')
            model_count = len(pdbid_chain_pdf['model'].unique())
            ca_count = pdbid_chain_pdf.shape[0]
            stats.append({
                'PDBid': pdbid,
                'het_hom': het_hom,
                'year': year,
                '#model': model_count,
                '#chains': total_chain_count,
                '#protChains': len(pdbid_chains_dict[pdbid]),
                'chain': chain,
                '#alphaCarbs': ca_count,
                'identity': identity,
                # 'similarity': 100  # TODO
            })
    pdf = pd.DataFrame(stats)



    pdf_sorted = pdf.sort_values(by=['year', '#model'], ascending=[True, True])
    print(pdf_sorted.dtypes)
    for col_to_cast in ['PDBid', 'chain']:
        pdf_sorted[col_to_cast] = pdf_sorted[col_to_cast].astype('string')
    print(pdf_sorted.dtypes)
    print(f'Completed {len(pdbid_chains)} CIFs in {round(time() - start)} seconds')
    return pdf_sorted


if __name__ == '__main__':
    # rpath_pdbid_chains = (f'../data/NMR/multimodel_PDBids/het_multimod_2104_pdbid_chains.txt')
    stats_pdf = generate_stats()

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

    # Note:
    # cif_dict = MMCIF2Dict.MMCIF2Dict(rpath_cif)
    # A few of these return '?', so using 'get_structure()' instead,
    # which broadly same but sometimes differs by as much as 3 years.
    # year = cif_dict.get('_citation.year', ['NA'])[0]
    # if year != 'NA' and year != '?':
    #     print(f'year is given as {year}. Replacing with 0')
    #     year = int(year)
    # else:
    #     year = 0

    # from Bio.PDB.Polypeptide import is_aa
    # first_model = next(structure.get_models())
    # for chain in first_model:
    #     ca_count = sum(1 for atom in chain.get_atoms() if atom.get_name() == 'CA')
    # num_models = len(list(structure))