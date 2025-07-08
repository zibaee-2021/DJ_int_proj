#!/usr/bin/env python3
"""
NOTE: THIS SCRIPT IS FROM https://www.dsimb.inserm.fr/ATLAS/data/download/download_ATLAS.txt

Bulk download script for ATLAS, Chameleon, and DPF MD datasets.
Supports fetching any combination of analysis, protein, total, or metadata archives per PDB chain.
"""
import sys, argparse, logging, time, shutil, zipfile, subprocess
from pathlib import Path

API_BASE = 'https://www.dsimb.inserm.fr/ATLAS/api'
DEFAULT_DEST = 'ATLAS_downloads'
LOCAL_DEST = '../../data/ATLAS_downloads'
SUPPORTED_DATASETS = ['ATLAS', 'chameleon', 'DPF']
AVAILABLE_ARCHIVES = ['analysis', 'protein', 'total', 'metadata']
MARKER = '.complete'


def check_aria2c_installed() -> None:
    """
    Ensure aria2c is available in PATH; exit if not.
    """
    if shutil.which('aria2c') is None:
        logging.error('aria2c is not installed or not found in PATH. Please install aria2c to use this script.')
        sys.exit(1)


def download_with_aria2c(url: str, output_path: Path, retries: int = 5, segments: int = 4, verify_ssl: bool = True) -> None:
    """
    Download a file with aria2c, enabling resume and segmentation.
    """
    cmd = [
        'aria2c',
        '--continue=true',
        f'--max-tries={retries}',
        f'--max-connection-per-server={segments}',
        '--remote-time=true',
        '--content-disposition=false',
        '-d', str(output_path.parent),
        '-o', output_path.name,
        url
    ]
    if not verify_ssl:
        cmd.extend(['--check-certificate=false'])
    logging.info(f'Running: {' '.join(cmd)}')
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f'aria2c failed with code {result.returncode}')


def download_parsable_lists(out_dir: Path, retries: int, verify_ssl: bool) -> dict[str, Path]:
    """
    Download unified parsable archive using aria2c, extract *_pdb.txt files,
    return mapping dataset->PDB list path.
    """
    url = f'{API_BASE}/parsable'
    logging.info(f'Fetching parsable lists from {url}')

    out_dir.mkdir(parents=True, exist_ok=True)
    final_zip = out_dir / 'parsable_latest.zip'

    # download with aria2c (resumable)
    download_with_aria2c(url, final_zip, retries=retries, verify_ssl=verify_ssl)

    mapping: dict[str, Path] = {}
    with zipfile.ZipFile(final_zip, 'r') as z:
        for name in z.namelist():
            if name.endswith('_pdb.txt'):
                ds = name.split('_')[-2]
                if ds in SUPPORTED_DATASETS:
                    path = out_dir / name
                    z.extract(name, path=out_dir)
                    mapping[ds] = path
    final_zip.unlink()

    if not mapping:
        logging.error('No PDB lists found for supported datasets.')
        sys.exit(1)
    return mapping


def read_pdb_list(txt_file: Path) -> list[str]:
    """
    Read non-empty chain IDs from a text file.
    """
    with open(txt_file) as f:
        ids = [line.strip() for line in f if line.strip()]
    logging.info(f'Loaded {len(ids)} identifiers from {txt_file}')
    return ids


def clean_incomplete(out_folder: Path):
    """
    Remove incomplete folder if marker not present.
    """
    marker = out_folder / MARKER
    if out_folder.exists() and not marker.exists():
        shutil.rmtree(out_folder)


def download_archive(dataset: str, pdb_chain: str, archive_type: str, dest: Path, retries: int, verify_ssl: bool) -> bool:
    """
    Download and extract a specific archive type for a PDB chain using aria2c.
    Handles metadata as JSON rather than ZIP.
    Returns True on success, False on failure.
    """
    out_folder = dest / dataset / pdb_chain / archive_type
    marker_file = out_folder / MARKER

    # Skip if already complete
    if marker_file.exists():
        return True

    # Clean any incomplete previous run
    clean_incomplete(out_folder)
    out_folder.mkdir(parents=True, exist_ok=True)

    endpoint = f'{API_BASE}/{dataset}/{archive_type}/{pdb_chain}'

    try:
        if archive_type == 'metadata':
            # Metadata is a JSON file, not a ZIP archive
            metadata_file = out_folder / f'{pdb_chain}_metadata.json'
            download_with_aria2c(endpoint, metadata_file, retries=retries, verify_ssl=verify_ssl)
            marker_file.write_text('done')
            return True
        else:
            # Other archives are ZIP files
            final_zip = dest / f'{dataset}_{pdb_chain}_{archive_type}.zip'
            # Download ZIP
            download_with_aria2c(endpoint, final_zip, retries=retries, verify_ssl=verify_ssl)
            # Extract
            with zipfile.ZipFile(final_zip, 'r') as z:
                z.extractall(out_folder)
            marker_file.write_text('done')
            final_zip.unlink()
            return True
    except Exception as e:
        logging.error(f'Failed to download {endpoint}: {e}')
        # Cleanup on failure
        if archive_type == 'metadata':
            if metadata_file.exists():
                metadata_file.unlink()
        else:
            if final_zip.exists():
                final_zip.unlink()
        if out_folder.exists():
            shutil.rmtree(out_folder)
        return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Bulk MD downloader with aria2c for resumable downloads.'
    )
    parser.add_argument(
        '--datasets', nargs='+', choices=SUPPORTED_DATASETS + ['all'],
        default=['ATLAS'],
        help='Datasets: ATLAS, chameleon, DPF, or all'
    )
    parser.add_argument(
        '--archives', nargs='+', choices=AVAILABLE_ARCHIVES + ['all'],
        default=['protein'],
        help='Archives: analysis, protein, total, metadata, or all'
    )
    parser.add_argument(
        '--pdb-list', type=Path,
        help='Text file of PDB IDs (e.g., 1abc_A) to filter download subset'
    )
    parser.add_argument(
        # '--dest', type=Path, default=Path(DEFAULT_DEST),
        '--dest', type=Path, default=Path(LOCAL_DEST),
        help='Root output directory'
    )
    parser.add_argument(
        '--retries', type=int, default=5,
        help='Retries on network errors'
    )
    parser.add_argument(
        '--segments', type=int, default=4,
        help='Number of parallel segments for aria2c'
    )
    parser.add_argument(
        '--verify-ssl', action='store_true',
        help='Enable SSL verification'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    # Ensure aria2c is available before proceeding
    check_aria2c_installed()

    dest = args.dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)

    mapping = download_parsable_lists(
        dest,
        retries=args.retries,
        verify_ssl=args.verify_ssl
    )
    user_chains = set(read_pdb_list(args.pdb_list)) if args.pdb_list else set()

    datasets = SUPPORTED_DATASETS if 'all' in args.datasets else args.datasets
    archives = AVAILABLE_ARCHIVES if 'all' in args.archives else args.archives

    # Build filtered chain lists and compute totals
    filtered: dict[str, list[str]] = {}
    for ds in datasets:
        if ds not in mapping:
            continue
        all_chains = read_pdb_list(mapping[ds])
        chains = [c for c in all_chains if c in user_chains] if user_chains else all_chains
        filtered[ds] = chains

    total_tasks = sum(len(chains) * len(archives) for chains in filtered.values())
    completed = 0
    stats = {ds: {arc: {'success':0,'fail':0} for arc in archives} for ds in filtered}

    for ds, chains in filtered.items():
        logging.info(f'Dataset {ds}: {len(chains)} chains, {len(archives)} archives each')
        for chain in chains:
            for arc in archives:
                ok = download_archive(
                    ds, chain, arc, dest,
                    retries=args.retries,
                    verify_ssl=args.verify_ssl
                )
                stats[ds][arc]['success' if ok else 'fail'] += 1
                completed += 1
                logging.info(f'Progress: {completed}/{total_tasks} ({ds}/{chain}/{arc})')
                time.sleep(0.2)

    # Final summary
    print('\nDownload Summary:')
    for ds, arcs in stats.items():
        print(f'Dataset: {ds}')
        for arc, counts in arcs.items():
            print(f'  {arc}: {counts['success']} succeeded, {counts['fail']} failed')
    print(f'\nCompleted {completed}/{total_tasks} tasks. Done.')


if __name__ == '__main__':
    main()
