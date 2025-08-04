import os, glob
from collections import defaultdict


def write_pdb_file_per_chain(rp_pdb_f: str, rp_dst_dir: str):
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

    with open(rp_pdb_f, 'r') as f:
        lines = f.readlines()
        # is there a pyhthon function that will return the index podition of a given string, e.g.  'MODEL        1'

    with open(rp_pdb_f, 'r') as f:
        for i, line in enumerate(f):
            record = line[:6].strip()
            if record in ['ATOM', 'TER', 'MODEL', 'ENDMDL']:
                if record == 'ATOM':
                    chain = line[21]
                    got_chain = True
                if got_chain and previous_record.startswith('MODEL'):
                    if previous_record != 'MODEL        1':
                        chain_lines[chain].append(f'{previous_record}\n')
                if got_chain and record.startswith('ATOM') and previous_record.startswith('TER'):  # 1 model can >1 chain
                    if chain != previous_chain:
                        chain_lines[previous_chain].append(f'ENDMDL\n')
                if got_chain:
                    chain_lines[chain].append(line)
                if line[:20].strip() == 'MODEL        1':
                    preamble = lines[0:i+1]
            previous_record = line[:20].strip()
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
    rp_pdb_dir = os.path.join('..', '..', 'data', 'NMR', 'raw_pdbs', 'hethom_combined')
    rp_pdbchains_dir = os.path.join('..', '..', 'data', 'NMR', 'pdb_chains', 'hethom_combined')
    os.makedirs(rp_pdbchains_dir, exist_ok=True)

    rp_pdb_files = sorted(glob.glob(os.path.join(rp_pdb_dir, '*.pdb')))

    for rp_pdb_f in rp_pdb_files:
        write_pdb_file_per_chain(rp_pdb_f, rp_pdbchains_dir)