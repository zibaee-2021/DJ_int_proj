import os, subprocess
import shutil
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def relpath_mmseqs_dir(het_hom: str) -> str:
    rp_mmseqs_dir = os.path.join('..', 'data', 'NMR', 'mmseqs', het_hom)
    return rp_mmseqs_dir


def write_fasta(het_hom: str, pdbid: str, pdbid_chains_dict: dict, bio_struct) -> str:
    rp_mmseqs_dir = relpath_mmseqs_dir(het_hom=het_hom)
    rp_fasta_dir = os.path.join(rp_mmseqs_dir, 'fasta')
    os.makedirs(rp_fasta_dir, exist_ok=True)

    model = bio_struct[0] # Using first model only
    ppb = PPBuilder()
    chain_sequences, records = {}, []
    for chain in model:
        if chain.id not in pdbid_chains_dict[pdbid]:
            continue
        else:
            peptides = ppb.build_peptides(chain)
            seq = ''.join(str(pp.get_sequence()) for pp in peptides)
            chain_sequences[chain.id] = seq
            record = SeqRecord(Seq(seq), id=f'{pdbid}_{chain.id}', description='')
            records.append(record)
    rp_fasta_file = os.path.join(rp_fasta_dir, f'{pdbid}_seqs.fasta')
    SeqIO.write(records, rp_fasta_file, 'fasta')
    return rp_fasta_file


def run_mmseqs_all_vs_all(rp_fasta_f, pdbid: str):
    """
    Runs MMseqs2 all-vs-all search on the provided FASTA file.
    Outputs:
      - results.m8 (tabular results)
      - Pandas DataFrame of alignments
    """
    het_hom = rp_fasta_f.split('/')[4]
    rp_mmseqs_dir = relpath_mmseqs_dir(het_hom=het_hom)
    rp_dir2delete = os.path.join(rp_mmseqs_dir, 'dir2delete')
    shutil.rmtree(rp_dir2delete, ignore_errors=True)

    os.makedirs(rp_dir2delete, exist_ok=True)
    tmp_dir = os.path.join(rp_dir2delete, f'{pdbid}_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    db_name = os.path.join(rp_dir2delete, f'{pdbid}_db')
    result_db = os.path.join(rp_dir2delete, f'{pdbid}_result')
    output_m8 = os.path.join(rp_dir2delete, f'{pdbid}_results.m8')

    # 1. Create database:
    subprocess.run(['mmseqs', 'createdb',
                    rp_fasta_f, db_name],
                   stdout=subprocess.DEVNULL, check=True)

    # 2. Search: all-vs-all:
    # subprocess.run(['mmseqs', 'search',
    #     db_name,   # query db
    #     db_name,   # target db
    #     result_db, tmp_dir, '--threads', '4'],
    #                stdout=subprocess.DEVNULL, check=True)

    # 2. Search: all-vs-all (with lower threshold):
    subprocess.run(['mmseqs', 'search',
                    db_name, db_name,
                    result_db, tmp_dir,
                    '--alignment-mode', '3',
                    '--threads', '4',
                    '-s', '7.5',
                    '-e', '1000'],
                   stdout=subprocess.DEVNULL, check=True)

    # 3. Convert alignments to tabular output
    subprocess.run(['mmseqs', 'convertalis',
                    db_name, db_name,
                    result_db, output_m8,
                    '--format-output',
                    'query,target,evalue,pident,alnlen'],
                   stdout=subprocess.DEVNULL, check=True)

    # 4. Read results into Pandas
    df = pd.read_csv(output_m8, sep='\t', header=None, names=['query','target','evalue','pident','alnlen'])
    return df


def filter_results(pdbid: str, het_hom: str, pdf):
    pdf = pdf[pdf['query'] != pdf['target']]

    rp_mmseqs_dir = relpath_mmseqs_dir(het_hom)

    if not pdf.empty:
        pdf_csv_results_dir = os.path.join(rp_mmseqs_dir, 'results')
        os.makedirs(pdf_csv_results_dir, exist_ok=True)
        pdf.to_csv(os.path.join(pdf_csv_results_dir, f'{pdbid}.csv'), index=False)
    else:
        print(f'No results for {pdbid}. Adding id to list file.')
        rp_zero_idty_lst = os.path.join(rp_mmseqs_dir, 'PDBid_no_idty.lst')

        if os.path.exists(rp_zero_idty_lst):
            with open(rp_zero_idty_lst, 'r') as f:
                pdbids_no_idty = f.readlines()
            pdbids_no_idty = [_pdbid.removesuffix('\n') for _pdbid in pdbids_no_idty]
            pdbids_no_idty.append(pdbid)
            pdbids_no_idty.sort()
            with open(rp_zero_idty_lst, 'w') as f:
                f.write('\n'.join(pdbids_no_idty) + '\n')
        else:
            with open(rp_zero_idty_lst, 'w') as f:
                f.write(f'{pdbid}\n')
    return pdf


if __name__ == '__main__':
    het_hom = 'heteromeric'
    # pdbid = '1A0N'
    # pdbid_chains_dict ={'1A0N': ['A', 'B']}
    pdbid = '1AOU'
    pdbid_chains_dict = {'1AOU': ['A', 'B']}
    rp_cif = os.path.join('..', 'data', 'NMR', 'raw_cifs', het_hom, f'{pdbid}.cif')
    parser = MMCIFParser(QUIET=True)
    bio_struct = parser.get_structure('', rp_cif)
    rp_fasta_file = write_fasta(het_hom, pdbid, pdbid_chains_dict, bio_struct)
    df_results = run_mmseqs_all_vs_all(rp_fasta_f=rp_fasta_file, pdbid=pdbid)
    print(df_results)
    filter_results(pdbid, het_hom, df_results)
    # df_results = df_results[df_results['query'] != df_results['target']]
    #
    # _rp_mmseqs_dir = relpath_mmseqs_dir(het_hom)
    # if not df_results.empty:
    #     _pdf_csv_results_dir = os.path.join(_rp_mmseqs_dir, 'results')
    #     os.makedirs(_pdf_csv_results_dir, exist_ok=True)
    #     df_results.to_csv(os.path.join(_pdf_csv_results_dir, f'{pdbid}.csv'), index=False)
    # else:
    #     print(f'No results for {pdbid}. Adding id to list file.')
    #     _rp_zero_idty_lst = os.path.join(_rp_mmseqs_dir, 'PDBid_no_idty.lst')
    #
    #     if os.path.exists(_rp_zero_idty_lst):
    #         with open(_rp_zero_idty_lst, 'r') as f:
    #             _pdbids_no_idty = f.readlines()
    #         _pdbids_no_idty = [_pdbid.removesuffix('\n') for _pdbid in _pdbids_no_idty]
    #         _pdbids_no_idty.append(pdbid)
    #         _pdbids_no_idty.sort()
    #         with open(_rp_zero_idty_lst, 'w') as f:
    #             f.write('\n'.join(_pdbids_no_idty) + '\n')
    #     else:
    #         with open(_rp_zero_idty_lst, 'w') as f:
    #             f.write(f'{pdbid}\n')

    pass