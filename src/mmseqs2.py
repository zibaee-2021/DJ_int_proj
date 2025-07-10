import os, glob
from time import time
import pandas as pd
from Bio.PDB import MMCIFParser, MMCIF2Dict, PDBParser, PPBuilder

cif_parser = MMCIFParser(QUIET=True)
pdb_parser = PDBParser(QUIET=True)


def _extract_aa_seq_from_pdb_or_cif(fpath, chain: str, cif_or_pdb: str) -> str:
    if cif_or_pdb == 'pdb':
        struct = pdb_parser.get_structure('', fpath)
    else:
        struct = cif_parser.get_structure('', fpath)
    chain_model = struct[0][chain]
    ppb = PPBuilder()
    aa_seq = ''
    for pp in ppb.build_peptides(chain_model):
        aa_seq = pp.get_sequence()
        print(aa_seq)
    return aa_seq


def _write_aa_seqs_to_fasta(pdbid_chain1: str, aa_seq1: str, pdbid_chain2: str, aa_seq2: str) -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    records = [SeqRecord(Seq(aa_seq1), id=pdbid_chain1), SeqRecord(Seq(aa_seq2), id=pdbid_chain2)]
    SeqIO.write(records, f'.fasta', 'fasta')


def calc_mmseq(fpath1, fpath2, cif_or_pdb: str):
    seq1 = _extract_aa_seq_from_pdb_or_cif(fpath1, cif_or_pdb)
    seq2 = _extract_aa_seq_from_pdb_or_cif(fpath2, cif_or_pdb)
    pdbid_chain1 = fpath1
    pdbid_chain2 = fpath2
    _write_aa_seqs_to_fasta(pdbid_chain1, seq1, pdbid_chain2, seq2)
    import subprocess
    # Create database (input_db, input_db.index, input_db.lookup, etc)
    subprocess.run(['mmseqs', 'createdb', f'{pdbid_chain1}_{pdbid_chain2}.fasta', 'input_db'], check=True)
    subprocess.run(['mmseqs', 'search',
                    'input_db', 'target_db', 'result_db',
                    'tmp_dir', '--threads', '4'], check=True)

    os.makedirs('tmp_dir', exist_ok=True)
    subprocess.run(['mmseqs', 'convertalis',
                    'input_db', 'target_db', 'result_db',
                    'results.m8', '--format-output', 'query,target,evalue,pident,alnlen'], check=True)

    df = pd.read_csv('results.m8', sep='\t', header=None, names=['query', 'target', 'evalue', 'pident', 'alnlen'])
    return df


if __name__ == '__main__':
    start = time()
    fpath1 = os.path.join('..', 'data', 'NMR', 'raw_cifs', 'heteromeric', '1A0N_A' )
    fpath2 = os.path.join('..', 'data', 'NMR', 'raw_cifs', 'heteromeric', '1A0N_B')
    pdf = calc_mmseq(fpath1, fpath2, cif_or_pdb='pdb')
    print(f'Completed ... in {round(time() - start)} seconds')