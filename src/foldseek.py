"""
Installation command for MacOS from FoldSeek github README.md:
`wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz`

Installed to my '~/bi_tools' directory and added `export PATH="$HOME/bi_tools/foldseek/bin:$PATH"` to my `.zprofile`
on 26Sep25. The version was: 8979d230fb64c7089380b652758d8705493ed4a5 (FoldSeek git commit hash hexidecimal SHA-1 digest).

Bug and workaround:
Attempting to run `foldseek convert2fasta {rp_prefix}_ss {rp_prefix}.3di.fasta` fails with missing header error.
This is a bug/issue reported in March 2022 (github.com/steineggerlab/foldseek/issues/15).
Workarounds provided at the time appear to still be needed (at least for my Mac version).

4-line workaround, provided by a FoldSeek contributor 'mvankem':

    mv db tmp
    cp db_ss db
    foldseek convert2fasta db db_ss.fasta
    mv tmp db

Alternative, and "better", 2-line workaround, provided by Steinegger lab member 'milot-mirdita':

    foldseek lndb queryDB_h queryDB_ss_h
    foldseek convert2fasta queryDB_ss queryDB_ss.fasta

"""
import subprocess
import os


def createdb4pdb(rp_prefix: str):
    """
    Run FoldSeek `createdb` function, writing 14 database files required for other downstream FoldSeek functionality.
    (Some files in binary format, so cannot be opened to view.)
    :param rp_prefix: Name of input pdb file and name/prefix of 14 output files.
    Must include corresponding relative path of input pdb file, which also serves as the destination path for all 14
    output files (e.g. '../data/XRAY/DynDom/Hayward_files/foldseek/').

    The 14 output database files:

        <prefix>: Binary database file containing the amino acid sequences (body).
        <prefix>.dbtype: Small file telling FoldSeek what type of file this database is.
        <prefix>.index: Index file mapping internal IDs to byte offsets in <prefix>, so entries can be fetched quickly.
        <prefix>.lookup: Maps internal numeric IDs back to the original identifiers (PDB chain IDs, etc.).
        <prefix>.source: Keeps track of where each entry came from (the input file paths, e.g. nameofpdb.pdb.

        <prefix>_ca: Database of CA-only coordinates extracted from structure. Used for fast structural prefiltering.
        <prefix>_ca.dbtype: Type marker for the CA database.
        <prefix>_ca.index: Index mapping IDs to the positions of the coordinate records inside <prefix>_ca

        <prefix>_h: Database containing headers/identifiers (what would appear after > in a FASTA,
                    i.e. separate from amino acid sequence 'body'.)
        <prefix>_h.dbtype: Type marker: tells Foldseek this DB holds headers, not bodies.
        <prefix>_h.index: Index file for fast lookups in the header database.

        <prefix>_ss: Database containing 3Di alphabet sequence for each residue; symbolic encoding of local 3D geometry.
        <prefix>_ss.dbtype: Type marker for the 3Di database.
        <prefix>_ss.index: Index for fast access to the 3Di database.
    """
    subprocess.run(['foldseek', 'createdb', f'{rp_prefix}.pdb', rp_prefix], check=True)

def _workaround(rp_prefix: str):
    """
    This workaround needs to be run before running `convert2fasta` for 3Di sequence, otherwise it will fail with
    'missing header' error message.
    It was provided by 'milot-mirdita' (Steinegger lab), see 'github.com/steineggerlab/foldseek/issues/15'.
    """
    subprocess.run(['foldseek', 'lndb', f'{rp_prefix}_h', f'{rp_prefix}_ss_h'], check=True)

def convert_3di_states_to_fasta(rp_prefix: str):
    """
    Write FoldSeek 3Di sequence to human-readable FASTA format (`.fasta`).
    This requires the workaround to be run first.
    :param rp_prefix:  Prefix of input and output filenames. (E.g. PDB filename).
    Must include corresponding relative path (e.g. '../data/XRAY/DynDom/Hayward_files/foldseek/'.), which is the same
    for both the input ('<rp_prefix>_ss') and output ('<rp_prefix>.3di.fasta') files.
    """
    _workaround(rp_prefix)
    subprocess.run(['foldseek', 'convert2fasta', f'{rp_prefix}_ss', f'{rp_prefix}.3di.fasta'], check=True)

def convert_aa_seq_to_fasta(rp_prefix: str):
    """
    Write amino acid sequence to FASTA format (`.fasta`).
    (This does NOT require any workaround to be run first.)
    :param rp_prefix: Name of input file and prefix of output filename. (E.g. PDB filename).
    Must include corresponding relative path (e.g. '../data/XRAY/DynDom/Hayward_files/foldseek/'.), which is the same
    for both the input ('<rp_prefix>') and output ('<rp_prefix>.aa.fasta') files.
    """
    subprocess.run(['foldseek', 'convert2fasta', rp_prefix, f'{rp_prefix}.aa.fasta'], check=True)


if  __name__ == '__main__':
    foldseek_dst_db_dir = '../data/DynDom/Hayward_files/foldseek/debug'
    os.makedirs(foldseek_dst_db_dir, exist_ok=True)
    rp_pdbname = os.path.join(foldseek_dst_db_dir, 'AACII1R1_confA_2a4n')
    # rp_pdbname = os.path.join(foldseek_dst_db_dir, 'AACII1R1_confB_1b87')
    createdb4pdb(rp_pdbname)  # generates 14 files
    # convert_3di_states_to_fasta(rp_pdbname)
    convert_aa_seq_to_fasta(rp_pdbname)
