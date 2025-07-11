import os, subprocess
import shutil
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def write_pdb_per_chain_to_fasta(fp_pdb: str, is_cif, rp_fasta):
    if is_cif:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure('structure', fp_pdb)
    model = structure[0] # I use model 0
    ppb = PPBuilder() # Init polypeptide builder
    chain_sequences = {}
    records = []

    for chain in model:  # Loop over all chains
        chain_id = chain.id

        # Some chains may have missing residues
        peptides = ppb.build_peptides(chain)
        if not peptides:
            continue

        # Concatenate all fragments for this chain
        seq = ''.join(str(pp.get_sequence()) for pp in peptides)
        chain_sequences[chain_id] = seq
        # Create SeqRecord:
        record = SeqRecord(
            Seq(seq),
            id=f'{os.path.basename(fp_pdb)}_{chain_id}',
            description=''
        )
        records.append(record)

    SeqIO.write(records, rp_fasta, 'fasta')  # Write FASTA
    return chain_sequences


def run_mmseqs_all_vs_all(rp_fasta, output_prefix: str, output_dir2delete: str):
    """
    Runs MMseqs2 all-vs-all search on the provided FASTA file.
    Outputs:
      - results.m8 (tabular results)
      - Pandas DataFrame of alignments
    """
    shutil.rmtree(output_dir2delete, ignore_errors=True)
    # Working subdirectories and files
    tmp_dir = os.path.join(output_dir2delete, f'{output_prefix}_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    db_name = os.path.join(output_dir2delete, f'{output_prefix}_db')
    result_db = os.path.join(output_dir2delete, f'{output_prefix}_result')
    output_m8 = os.path.join(output_dir2delete, f'{output_prefix}_results.m8')

    # 1. Create database:
    subprocess.run(['mmseqs', 'createdb',
                    rp_fasta, db_name],
                   stdout=subprocess.DEVNULL, check=True)

    # 2. Search: all-vs-all:
    # subprocess.run(['mmseqs', 'search',
    #     db_name,   # query db
    #     db_name,   # target db
    #     result_db, tmp_dir, '--threads', '4'],
    #                stdout=subprocess.DEVNULL, check=True)

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


if __name__ == '__main__':
    pdbid = '1A0N'
    _meric = 'heteromeric'
    relpath_pdb = f'../data/NMR/raw_cifs/{_meric}/{pdbid}.cif'  # or '/1A0N.pdb'
    mmseqs_dir = '../data/NMR/mmseqs'
    fasta_dir = os.path.join(mmseqs_dir, 'fasta')
    os.makedirs(fasta_dir, exist_ok=True)
    relpath_fasta = os.path.join(fasta_dir, f'{pdbid}_seqs.fasta')
    seqs = write_pdb_per_chain_to_fasta(fp_pdb=relpath_pdb, is_cif=relpath_pdb.endswith('.cif'), rp_fasta=relpath_fasta)
    print('Extracted sequences:')
    for chain, sequence in seqs.items():
        print(f'Chain {chain}: {sequence}')

    output_dir2delete = os.path.join(mmseqs_dir, 'dir2delete')
    os.makedirs(output_dir2delete, exist_ok=True)
    output_dir2delete = os.path.join(output_dir2delete, pdbid)
    df_results = run_mmseqs_all_vs_all(rp_fasta=relpath_fasta, output_prefix=_meric, output_dir2delete=output_dir2delete)
    print(df_results)
    pdf_csv_results_dir = os.path.join(mmseqs_dir, 'results')
    os.makedirs(pdf_csv_results_dir, exist_ok=True)
    df_results.to_csv(os.path.join(pdf_csv_results_dir, f'{pdbid}.csv'), index=False)
    pass
