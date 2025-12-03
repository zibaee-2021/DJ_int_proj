"""
Inspired by:
'Rigid Domains in Proteins: An Algorithmic Approach to Their Identification'
Nichols W.L., Rose G.D, Eyk L.F.T., Zimm B.H.
Proteins: Struct., Func., and Gen. 23:38-48 (1995)
"""
import os, glob
from typing import List, Tuple
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score
from Bio.SVDSuperimposer import SVDSuperimposer

def _rp_nmr_dir() -> str:
    return os.path.join('..', 'data', 'NMR')

def _rp_parsed_cifs(subdir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', subdir)

def _distance_matrix(coords: np.ndarray) -> np.ndarray:
    """
    N×N Euclidean distance matrix for one model.
    D[i,j] = ||coords[i] - coords[j]||_2

    To calculate Euclidean distance between each residue with every other in given 2D array of coords:
    - Convert (N x (x,y,z)) coordinates to column vector; and to row vector.
    - Subtract column vector by row vector, giving the distance matrix of each x-x, y-y, z-z.
    - Convert to Euclidean distances by taking norm of these per coord distances, giving an NxN distance matrix.
    """
    # Using None in this way is just another way (than e.g. using `.reshape()`) for reshaping coords (N,3) to diff (N,N,3).
    coords_col_vec = coords[:, None, :]  # (N,1,3) <- (N,3)
    coords_row_vec = coords[None, :, :]  # (1,N,3) <- (N,3)
    diff = coords_col_vec - coords_row_vec  # (N,N,3) <- (N,1,3) - (1,N,3)
    dist_matrx = np.linalg.norm(diff, axis=-1)  # (N,N) <- (N,N,3)
    np.fill_diagonal(dist_matrx, 0.0)  # should already have 0.0 along diagonal, but just reinforces this.
    return dist_matrx

def _align(coords_mod_i, coords_mod_j):
    assert len(coords_mod_i) == len(coords_mod_j), 'Sequences should be of identical length.'
    si = SVDSuperimposer()
    si.set(reference_coords=coords_mod_i, coords=coords_mod_j)
    si.run()
    return si

def compute_ddms(pidc_pdf):
    ddms = []
    model_numbers = pidc_pdf['A_pdbx_PDB_model_num'].unique()
    model_numbers = list(model_numbers)
    model_pairs = itertools.combinations(model_numbers, 2)
    mp_list = []
    for model1, model2 in model_pairs:
        mp_list.append((int(model1), int(model2)))
        model_i_pdf = pidc_pdf[pidc_pdf['A_pdbx_PDB_model_num'] == model1].copy()
        model_j_pdf = pidc_pdf[pidc_pdf['A_pdbx_PDB_model_num'] == model2].copy()
        coords_i = model_i_pdf[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
        coords_j = model_j_pdf[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values
        dm_i = _distance_matrix(coords_i)
        dm_j = _distance_matrix(coords_j)
        ddm = np.abs(dm_i - dm_j)
        np.fill_diagonal(ddm, 0.0)  #
        ddms.append(ddm)
    return ddms, mp_list


def compute_ddm_of_pair(pidc_pdf, pidc_pdf2):
    coords_i = pidc_pdf[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values  # 2DN1
    coords_j = pidc_pdf2[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']].values  # 2DN2
    coords_j = coords_j[1:][:]
    si = _align(coords_i, coords_j)
    model_i_si, model_j_si = si.reference_coords, si.coords
    dm_i = _distance_matrix(model_i_si)
    dm_j = _distance_matrix(model_j_si)
    ddm = np.abs(dm_i - dm_j)
    np.fill_diagonal(ddm, 0.0)  #
    return ddm

def analyse_ddms(ddms: list, model_pairs: List[tuple]=None, min_seq_sep: int=5, k: int=None, top_k_pairs: int=20,
                 random_state: int=0):
    """
    Analyse list of difference distance matrices (ddms) for a single PDB-chain.
    If k=None, auto-choose k in {2,3,4} by silhouette.

    :param ddms: Each element is an (N,N) symmetric matrix with zero diagonal: ddm = dm(model_b) - dm(model_a).
    :param model_pairs: Labels for each DDM (e.g., (model_i, model_j)). If None, indices 0..K-1 are used.
    :param min_seq_sep: Ignore near-diagonal |i-j| < min_seq_sep when computing scores (suppresses local noise).
    :param k: Number of residue domains. If None, choose best in {2,3,4} by silhouette on features.
    :param top_k_pairs: How many model pairs to return as "most different".
    :param random_state: For deterministic clustering.
    :return: results:
                    domain_labels : (N,) np.ndarray[int]
                    hinge_score   : (N,) np.ndarray[float]
                    per_res_max   : (N,) np.ndarray[float]
                    pair_scores   : list[(pair_label, float)], sorted desc
                    agg_rms       : (N,N) np.ndarray  # RMS across pairs of |ddm|
                    agg_max       : (N,N) np.ndarray  # max across pairs of |ddm|
                    similarity    : (N,N) np.ndarray  # residue×residue cosine similarity used for clustering
                    mask          : (N,N) boolean     # long-range mask used
    """
    if len(ddms) == 0:
        raise ValueError("ddms list is empty.")
    N = ddms[0].shape[0]
    for ddm in ddms:
        if ddm.shape != (N, N):
            raise ValueError("All DDMs must have identical (N,N) shape.")
        if not np.allclose(ddm, ddm.T, atol=1e-6):
            raise ValueError("DDM not symmetric within tolerance.")
        if not np.allclose(np.diag(ddm), 0.0, atol=1e-6):
            raise ValueError("DDM diagonal must be ~0.")

    # Long-range mask
    i, j = np.ogrid[:N, :N]
    mask = (np.abs(i - j) >= min_seq_sep)
    np.fill_diagonal(mask, False)

    # K = len(ddms)
    abs_ddm = np.stack([np.abs(ddm) for ddm in ddms], axis=0)  # (K, N, N)

    # Aggregate over all pairs
    agg_rms = np.sqrt(np.mean(abs_ddm**2, axis=0))         # RMS across pairs
    agg_max = np.max(abs_ddm, axis=0)                      # max across pairs

    # Build per-residue feature vectors from long-range parts of aggregates
    # Feature = [ agg_rms[i, long-range j]; agg_max[i, long-range j] ]
    # Zero out near-diagonal (short-range) entries but keep full length per row
    agg_rms_masked = agg_rms * mask  # (N,N)
    agg_max_masked = agg_max * mask  # (N,N)
    # Feature for residue r is the concatenation of the full masked rows (fixed length 2N)
    X = np.hstack([agg_rms_masked, agg_max_masked])  # (N,2N)

    # Residue similarity (cosine on features)
    X_abs = np.abs(X)  # (N,2N)                   # ignore sign; we're about magnitude patterns
    similarity = cosine_similarity(X_abs)  # (N,N)
    np.fill_diagonal(similarity, 1.0)
    similarity[similarity < 0] = 0.0  # clamp tiny negatives

    # Choose k if needed (by silhouette on X)
    if k is None:
        best_k, best_labels, best_s = None, None, -np.inf
        for kk in (2, 3, 4):
            labels_try = SpectralClustering(
                n_clusters=kk, affinity='precomputed',
                assign_labels='kmeans', random_state=random_state
            ).fit_predict(similarity)
            # Silhouette on the original feature space X (Euclidean)
            # If a cluster collapses (rare), silhouette can fail; guard it.
            try:
                s = silhouette_score(X, labels_try, metric='euclidean')
            except Exception:
                s = -1e9
            if s > best_s:
                best_k, best_labels, best_s = kk, labels_try, s
        k = best_k
        domain_labels = best_labels
    else:
        domain_labels = SpectralClustering(
            n_clusters=k, affinity='precomputed',
            assign_labels='kmeans', random_state=random_state
        ).fit_predict(similarity)

    # Per-residue "who moved?" summaries for visualisation
    # 1) max over partners and pairs (strong, localised differences)
    per_residue_max = np.max(agg_max * mask, axis=1)
    # 2) overall RMS over partners and pairs (global participation in change)
    per_residue_rms = np.sqrt(np.mean((agg_rms**2) * mask, axis=1))

    # Hinge score = z(per_residue_rms) * boundaryness
    # boundaryness: sequential change of domain label (smoothed)
    boundary = np.zeros(N)
    boundary[1:] += (domain_labels[1:] != domain_labels[:-1]).astype(float)
    boundary[:-1] += (domain_labels[:-1] != domain_labels[1:]).astype(float)
    boundary *= 0.5
    # light smoothing
    kernel = np.array([0.25, 0.5, 0.25])
    boundary_smooth = np.convolve(boundary, kernel, mode='same')

    # z-score of per_residue_rms (robust)
    med = np.median(per_residue_rms)
    mad = np.median(np.abs(per_residue_rms - med)) + 1e-12
    z_rms = (per_residue_rms - med) / (1.4826 * mad)
    hinge_score = (z_rms.clip(min=0) + 1.0) * boundary_smooth  # emphasise boundaries with high motion

    # Rank model pairs by global difference (RMS over long-range entries)
    denom = np.sqrt(mask.sum())
    pair_scores = []
    for idx, ddm in enumerate(ddms):
        score = np.linalg.norm(ddm[mask]) / denom
        label = model_pairs[idx] if model_pairs is not None else idx
        pair_scores.append((label, float(score)))
    pair_scores.sort(key=lambda t: t[1], reverse=True)
    if top_k_pairs is not None:
        pair_scores = pair_scores[:top_k_pairs]

    return dict(
        domain_labels=domain_labels,
        hinge_score=hinge_score,
        per_residue_max=per_residue_max,
        pair_scores=pair_scores,
        agg_rms=agg_rms,
        agg_max=agg_max,
        similarity=similarity,
        mask=mask,
        k=k,
    )

def analyse_single_ddm(ddm, model_pairs: List[tuple]=None, min_seq_sep: int=5, k: int=None, top_k_pairs: int=20,
                       random_state: int=0):
    ddms = [ddm]
    N = ddms[0].shape[0]
    if ddm.shape != (N, N):
        raise ValueError("All DDMs must have identical (N,N) shape.")
    if not np.allclose(ddm, ddm.T, atol=1e-6):
        raise ValueError("DDM not symmetric within tolerance.")
    if not np.allclose(np.diag(ddm), 0.0, atol=1e-6):
        raise ValueError("DDM diagonal must be ~0.")

    # Long-range mask
    i, j = np.ogrid[:N, :N]
    mask = (np.abs(i - j) >= min_seq_sep)
    np.fill_diagonal(mask, False)

    abs_ddm = np.stack([np.abs(ddm) for ddm in ddms], axis=0)  # (K, N, N)

    # Aggregate over all pairs
    agg_rms = np.sqrt(np.mean(abs_ddm**2, axis=0))         # RMS across pairs
    agg_max = np.max(abs_ddm, axis=0)                      # max across pairs

    # Build per-residue feature vectors from long-range parts of aggregates
    # Feature = [ agg_rms[i, long-range j]; agg_max[i, long-range j] ]
    # Zero out near-diagonal (short-range) entries but keep full length per row
    agg_rms_masked = agg_rms * mask  # (N,N)
    agg_max_masked = agg_max * mask  # (N,N)
    # Feature for residue r is the concatenation of the full masked rows (fixed length 2N)
    X = np.hstack([agg_rms_masked, agg_max_masked])  # (N,2N)

    # Residue similarity (cosine on features)
    X_abs = np.abs(X)  # (N,2N)                   # ignore sign; we're about magnitude patterns
    similarity = cosine_similarity(X_abs)  # (N,N)
    np.fill_diagonal(similarity, 1.0)
    similarity[similarity < 0] = 0.0  # clamp tiny negatives

    # Choose k if needed (by silhouette on X)
    if k is None:
        best_k, best_labels, best_s = None, None, -np.inf
        for kk in (2, 3, 4):
            labels_try = SpectralClustering(
                n_clusters=kk, affinity='precomputed',
                assign_labels='kmeans', random_state=random_state
            ).fit_predict(similarity)
            # Silhouette on the original feature space X (Euclidean)
            # If a cluster collapses (rare), silhouette can fail; guard it.
            try:
                s = silhouette_score(X, labels_try, metric='euclidean')
            except Exception:
                s = -1e9
            if s > best_s:
                best_k, best_labels, best_s = kk, labels_try, s
        k = best_k
        domain_labels = best_labels
    else:
        domain_labels = SpectralClustering(
            n_clusters=k, affinity='precomputed',
            assign_labels='kmeans', random_state=random_state
        ).fit_predict(similarity)

    # Per-residue "who moved?" summaries for visualisation
    # 1) max over partners and pairs (strong, localised differences)
    per_residue_max = np.max(agg_max * mask, axis=1)
    # 2) overall RMS over partners and pairs (global participation in change)
    per_residue_rms = np.sqrt(np.mean((agg_rms**2) * mask, axis=1))

    # Hinge score = z(per_residue_rms) * boundaryness
    # boundaryness: sequential change of domain label (smoothed)
    boundary = np.zeros(N)
    boundary[1:] += (domain_labels[1:] != domain_labels[:-1]).astype(float)
    boundary[:-1] += (domain_labels[:-1] != domain_labels[1:]).astype(float)
    boundary *= 0.5
    # light smoothing
    kernel = np.array([0.25, 0.5, 0.25])
    boundary_smooth = np.convolve(boundary, kernel, mode='same')

    # z-score of per_residue_rms (robust)
    med = np.median(per_residue_rms)
    mad = np.median(np.abs(per_residue_rms - med)) + 1e-12
    z_rms = (per_residue_rms - med) / (1.4826 * mad)
    hinge_score = (z_rms.clip(min=0) + 1.0) * boundary_smooth  # emphasise boundaries with high motion

    # Rank model pairs by global difference (RMS over long-range entries)
    denom = np.sqrt(mask.sum())
    pair_scores = []
    for idx, ddm in enumerate(ddms):
        score = np.linalg.norm(ddm[mask]) / denom
        label = model_pairs[idx] if model_pairs is not None else idx
        pair_scores.append((label, float(score)))
    pair_scores.sort(key=lambda t: t[1], reverse=True)
    if top_k_pairs is not None:
        pair_scores = pair_scores[:top_k_pairs]

    return dict(
        domain_labels=domain_labels,
        hinge_score=hinge_score,
        per_residue_max=per_residue_max,
        pair_scores=pair_scores,
        agg_rms=agg_rms,
        agg_max=agg_max,
        similarity=similarity,
        mask=mask,
        k=k,
    )

def plot_sequence_signal(signal, title='Per-residue signal', ylabel='magnitude', save_png=None):
    fig, ax = plt.subplots(figsize=(10, 2.8))
    ax.plot(np.arange(1, signal.shape[0] + 1), signal, linewidth=1.5)
    ax.set_xlabel('Residue index')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.margins(x=0)
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
    if save_png:
        plt.savefig(save_png, dpi=200, bbox_inches='tight')
    plt.show()
    return fig, ax

def plot_matrix(M, title='Matrix', save_png=None):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(M, origin='upper', aspect='auto')
    ax.set_title(title)
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Residue index')
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if save_png:
        plt.savefig(save_png, dpi=200, bbox_inches='tight')
    plt.show()
    return fig, ax


if __name__ == "__main__":

    look_at_haem_xray_pair_only = False
    if look_at_haem_xray_pair_only:
        # DDM BETWEEN HAEMOGLOBIN B-CHAIN OXY and DEOXY:
        rp_xray_dir = os.path.join('..', 'data', 'XRAY', 'handpicked')
        rp_haem_cifs = sorted(glob.glob(os.path.join(rp_xray_dir, 'haemoglobin', '**', '*.cif'), recursive=True))
        rp_pidc_ssvs = sorted(glob.glob(os.path.join(rp_xray_dir, 'parsed_cifs', '*.ssv'), recursive=True))
        pidc_pdf = pd.read_csv(rp_pidc_ssvs[0], sep=' ')
        pidc_pdf2 = pd.read_csv(rp_pidc_ssvs[1], sep=' ')
        _ddm = compute_ddm_of_pair(pidc_pdf, pidc_pdf2)
        res = analyse_single_ddm(_ddm, min_seq_sep=5)

        plot_sequence_signal(res['per_residue_max'], title='Max |ddm| per residue (any partner/pair)', ylabel='Å')
        plot_sequence_signal(res['hinge_score'], title='Hinge score per residue', ylabel='score')
        plot_matrix(res['agg_rms'], title='RMS( |ddm| ) across model pairs')

    rp_pidc_parsed_cif_dir = _rp_parsed_cifs('multimod_2713_hetallchains_hom1chain')
    rp_pidc_ssvs = sorted(glob.glob(os.path.join(rp_pidc_parsed_cif_dir, '*.ssv')))
    for rp_pidc_ssv in  rp_pidc_ssvs:
        pidc = os.path.basename(rp_pidc_ssv).removesuffix('.ssv')
        print(pidc)
        _pdf = pd.read_csv(rp_pidc_ssv, sep=' ')
        _ddms, _model_pairs = compute_ddms(_pdf)
        res = analyse_ddms(_ddms, model_pairs=_model_pairs, min_seq_sep=5)

        domain_labels = res['domain_labels']        # (N,)
        hinge_score = res['hinge_score']            # (N,)
        permax = res['per_residue_max']             # (N,)
        top_pairs = res['pair_scores']              # list of (pair, score) (PDB with 4 models -> 6 pairs -> list of 6).
        agg_rms = res['agg_rms']                    # (N,N)
        agg_max = res['agg_max']                    # (N,N)
        similarity = res['similarity']              # (N,N)
        mask = res['mask']                          # (N,N)
        k = res['k']                                # 2, 3, or 4

        # (a) Nice sequence visualisations:
        plot_sequence_signal(res['per_residue_max'], title='Max |ddm| per residue (any partner/pair)', ylabel='Å')
        plot_sequence_signal(res['hinge_score'], title='Hinge score per residue', ylabel='score')
        plot_matrix(res['agg_rms'], title='RMS( |ddm| ) across model pairs')

        # (b) Pairs that are likely different:
        top_scoring_pairs = []
        for pair, score in top_pairs:
            top_scoring_pairs.append(f'{pair} {score:.3f}')
            # print(f'{pair} {score:.3f}')
        rp_prot_dynamics_dir = os.path.join(_rp_nmr_dir(), 'prot_dynamics')
        os.makedirs(rp_prot_dynamics_dir, exist_ok=True)
        with open(os.path.join(rp_prot_dynamics_dir, f'{pidc}.lst'), 'w') as f:
            f.write('\n'.join(top_scoring_pairs))


