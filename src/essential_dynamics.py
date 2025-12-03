"""
Inspired by:
'Principal Component Analysis: A Method for Determining the Essential Dynamics of Proteins'
Charles C. David and Donald J. Jacobs
Methods Mol Biol. 1084: 193–226. (2014)
----------------------------------------------------------------------------------------------------------
Originates in 1980s onwards (in terms of proteins):
Method for Estimating the Configurational Entropy of Macromolecules.
Karplus M. & Kushick J.N.
Macromolecules 14: 325-332. (1981)

The effects of solvent on the conformation and the collective motions of protein: NMA and MD simulations of melittin in water and in vacuum.
Kitao, A., Hirata F., Go, N.
Chemical Physics 158: 447-472. (1991)

Large-amplitude nonlinear motions in proteins.
Garcia, A.E.
Phys. Rev. Lett. 68:2696–2699. (1992)

Effect of solvent on collective motions in globular proteins.
Hayward S., Kitao, A., Hirata, F., Go, N.
J. Mol. Biol. 234:1207–1217. (1993)
----------------------------------------------------------------------------------------------------------
"""

import itertools
import os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import KernelPCA

def _rp_nmr_dir() -> str:
    return os.path.join('..', 'data', 'NMR')

def _rp_parsed_cifs(subdir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', subdir)

def kabsch(P, Q):
    """Align P (N,3) to Q (N,3) using Kabsch; returns aligned P."""
    Pc = P - P.mean(axis=0, keepdims=True)
    Qc = Q - Q.mean(axis=0, keepdims=True)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    return Pc @ R + Q.mean(axis=0, keepdims=True)

def align_all_to_ref(coords, ref_idx=0):
    """
    Align all models in coords (M, N, 3) to coords[ref_idx] using Kabsch.
    """
    M, _, _ = coords.shape
    ref = coords[ref_idx]
    out = coords.copy()
    for i in range(M):
        out[i] = kabsch(coords[i], ref)
    return out

def build_data_matrix(aligned_coords, use_correlation=False):
    """Return (A_centered, mean_per_row, row_std) with A having shape (3N, M)."""
    M, N, _ = aligned_coords.shape
    A = aligned_coords.reshape(M, 3 * N).T  # (3N, M), columns are models (observations)
    row_means = A.mean(axis=1, keepdims=True)
    A0 = A - row_means
    if use_correlation:
        row_std = A0.std(axis=1, ddof=1, keepdims=True)
        row_std[row_std == 0] = 1.0
        A0 = A0 / row_std
    else:
        row_std = None
    return A0, row_means, row_std

def pca_evd(A0):
    """
    PCA via EVD of covariance Q = A0 @ A0.T;
    returns eigvals, eigvecs, scores.
    """
    # Q is (3N x 3N); With typical N this is fine. But for very large 3N, might be better to switch to SVD on (M x 3N).
    Q = A0 @ A0.T
    # Covariance matrix is always guaranteed to be symmetric, so always diagonalisable, so eigenvalues are guaraneteed
    # to all be real, and eigenvectors are guaranteed to all be orthogonal - which is needed for PCA.
    # We know Q is symmetric so we use np.linalg.eigh() instead of np.linalg.eig().
    eigvals, eigvecs = np.linalg.eigh(Q)

    # PCA requires ordered eigenvectors by decreasing variance (and eigendecomposition does not guarantee order), so:
    order = np.argsort(eigvals)[::-1]  # orders by smallest-to-largest), [::-1] reverses it to largest-to-smallest.
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    # PC scores for each model: project each centered column (i.e. model) onto each PCA mode (i.e. eigenvectors).
    scores = eigvecs.T @ A0  # shape: (3N_modes, M)
    return eigvals, eigvecs, scores  # scores[i] is PC_i across models

# NOTE: `pca_evd()` above is a manual implementation of PCA which has a few advantages over
# e.g. scikit-learn PCA which could be implemented as below (though not currently used):
from sklearn.decomposition import PCA
def pca_sklearn(A0, n_components=None):
    # A0 is (3N, M) with zero mean per row
    X = A0.T  # shape: (M, 3N), models × features
    pca = PCA(n_components=n_components, svd_solver='full')
    scores = pca.fit_transform(X)  # (M, n_components)
    eigvals = pca.explained_variance_  # eigenvalues
    eigvecs = pca.components_.T  # (3N, n_components)
    return eigvals, eigvecs, scores.T  # transpose scores for consistency

def weighted_rmsd_modes(eigvecs, eigvals, N, modes=(0, 1, 2)):
    """
    Map 3N eigenvector components -> per-residue RMSD contributions,
    weighted by sqrt(eigenvalue) to get Å units (Recipe I step 15).
    """
    out = {}
    for mi in modes:
        v = eigvecs[:, mi].reshape(N, 3)
        # per-residue amplitude (Å): sqrt( sum_{x,y,z} v^2 ) * sqrt(lambda)
        amp = np.linalg.norm(v, axis=1) * np.sqrt(max(eigvals[mi], 0))  #norm gives shape of motion per residue.
        out[mi] = amp
    return out  # dict: mode_index -> (N,)

def plot_scree(eigvals, title='Scree (eigenvalues)'):
    cum = np.cumsum(eigvals) / np.sum(eigvals)
    fig, ax1 = plt.subplots(figsize=(6.2, 3.6))
    ax1.plot(np.arange(1, len(eigvals) + 1), eigvals, marker='o', lw=1)
    ax1.set_xlabel('Mode index')
    ax1.set_ylabel('Eigenvalue (variance)')
    ax1.set_title(title)
    ax1.grid(True, ls='--', lw=0.5, alpha=0.5)
    ax2 = ax1.twinx()
    ax2.plot(np.arange(1, len(eigvals) + 1), cum, marker='.', lw=1)
    ax2.set_ylabel('Cumulative variance')
    plt.show()
    return fig, (ax1, ax2)

def scatter_pc(scores, i=0, j=1, labels=None, title=''):
    x, y = scores[i], scores[j]
    fig, ax = plt.subplots(figsize=(5, 4.2))
    if labels is None:
        ax.scatter(x, y, s=20, alpha=0.85)
    else:
        for lab in np.unique(labels):
            sel = (labels == lab)
            ax.scatter(x[sel], y[sel], s=26, alpha=0.9, label=str(lab))
        ax.legend(frameon=False)
    ax.set_xlabel(f'PC{i + 1}')
    ax.set_ylabel(f'PC{j + 1}')
    ax.set_title(title)
    ax.grid(True, ls='--', lw=0.5, alpha=0.5)
    plt.show()
    return fig, ax

# (Reuse your helpers for per-residue mode plots)
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

def essential_dynamics_pca(pidc: str, pidc_pdf, use_correlation=False, top_modes_for_rmsd=(0, 1, 2),
                           try_kmeans_k=None, pca_top_for_kpca=5, do_kpca=False, kpca_gamma=None):
    """
    pidc: PDB-chain identifier.
    pidc_pdf: PDB-chain dataframe.
    use_correlation: if True, normalize per-variable before PCA (paper advises trying both).
    top_modes_for_rmsd: tuple of mode indices for per-residue hinge profiles.
    try_kmeans_k: if int, run k-means on PC scores and print silhouette.
    do_kpca: if True, run KernelPCA (RBF) on the first pca_top_for_kpca PCs for extra separation.
    kpca_gamma: float ≥ 0 or None. If None (or if you previously used 'scale'/'auto'),
    we set gamma = 1.0 / n_features on the standardised top-PCs (a simple heuristic).
    """
    # Ensure sorting by model then row index within each model
    # df_sorted = pidc_pdf.sort_values(['A_pdbx_PDB_model_num']).reset_index(drop=True)

    # Group by model number and stack coordinates.
    # E.g. for 1A03_A which has 20 models, each with 90 CAs and 7 cols of data:
    # pidc_pdf has shape (1800,7) -> np array shape (20,) containing 20 arrays of shape (90,3) -> np array (20,90,3).
    # coords = (df_sorted
    coords = (pidc_pdf
              .groupby('A_pdbx_PDB_model_num')[['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']]
              .apply(lambda g: g.to_numpy())
              .to_numpy()
              )
    # coords are CA coordinates across models.
    # shape (M, N, 3) where M is number of models, N is number of atoms (CAs), 3 are x,y,z coordinates.
    # E.g. PDB with 20 models, 90 CAs has shape (20, 90, 3).
    coords = np.stack(coords, axis=0)

    M, N, _ = coords.shape
    aligned = align_all_to_ref(coords, ref_idx=0)

    # Build the data matrix for PCA from aligned coordinates.
    # `use_correlation=True` standardises per-coordinate (correlation PCA); otherwise covariance PCA.
    co_PCA, _, _ = build_data_matrix(aligned, use_correlation)

    # Compute eigendecomposition/SVD-based PCA of the covariance/correlation.
    # Returns eigenvalues (mode variances), eigenvectors (mode shapes), and scores (model projections).
    eigvals, eigvecs, scores = pca_evd(co_PCA)  # eigvals: (3N'), eigvecs: (3N', 3N'), scores: (3N', M)

    # For manual visual analysis:
    # Plot variance distribution across PCs to pick a meaningful cutoff (“elbow”/scree test).
    plot_scree(eigvals, title=f"{pidc} Scree ({'corr' if use_correlation else 'cov'}-PCA)")

    # For manual visual analysis:
    # Inspect separation of models in low-dimensional PC space (PC1 vs PC2) for clustering tendencies.
    scatter_pc(scores, 0, 1, title=f'{pidc} Models in PC1–PC2 space')

    # Optional quick unsupervised clustering directly in PC space to label putative ensembles.
    # We take the first few PCs, cluster models with k-means, compute silhouette as a separation quality metric,
    # and re-plot PC1–PC2 colored by cluster for interpretability.
    if try_kmeans_k:
        X = scores[:max(2, min(10, scores.shape[0])), :].T  # models × PCs
        km = KMeans(n_clusters=try_kmeans_k, n_init='auto', random_state=0).fit(X)
        km_labels = km.labels_
        s = silhouette_score(X, km_labels) if try_kmeans_k > 1 else np.nan
        print(f'[k-means] k={try_kmeans_k}, silhouette={s:.3f}')
        scatter_pc(scores, 0, 1, labels=km_labels, title=f'{pidc} PC1–PC2 with k-means (k={try_kmeans_k})')

    # Per-residue weighted RMSD modes (finds high amplitude residues; may be flexible regions (i.e. hinge/loop-like).
    n_used = N
    wr = weighted_rmsd_modes(eigvecs, eigvals, n_used, modes=top_modes_for_rmsd)
    for mi, amp in wr.items():
        title = f'{pidc} Weighted RMSD Mode (PC{mi + 1})'
        plot_sequence_signal(amp, title=title, ylabel='Å (weighted)')

    # Optional: kernel PCA on top of standard PCs (Recipe III)
    kpca_scores = None
    if do_kpca:
        # compress with standard PCA first (authors show using ~5 PCs works well)
        X = scores[:pca_top_for_kpca, :].T  # models × top PCs
        # standardise those PCs before RBF kernel for balance
        Xz = StandardScaler().fit_transform(X)

        # --- gamma handling ---
        # If user gives None (or used 'scale'/'auto' before), choose a simple heuristic
        # analogous to SVM's 'scale': gamma = 1 / n_features on standardised PCs—robust default without extra tuning.
        if kpca_gamma is None or (isinstance(kpca_gamma, str) and kpca_gamma.lower() in {'scale', 'auto'}):
            gamma = 1.0 / Xz.shape[1]
        elif isinstance(kpca_gamma, (int, float)) and kpca_gamma >= 0:
            gamma = float(kpca_gamma)
        else:
            raise ValueError("kpca_gamma must be None, a non-negative float, or 'scale'/'auto'.")

        # Compute a 2D non-linear embedding capturing curved manifolds of motion in model space.
        kpca = KernelPCA(n_components=2, kernel='rbf', gamma=gamma, fit_inverse_transform=False)
        Xk = kpca.fit_transform(Xz)
        kpca_scores = Xk.T

        # Visualise non-linear separation; clusters that overlap linearly may separate here.
        fig, ax = plt.subplots(figsize=(5, 4.2))
        ax.scatter(Xk[:, 0], Xk[:, 1], s=24, alpha=0.9)
        ax.set_xlabel('kPC1')
        ax.set_ylabel('kPC2')
        ax.set_title(f'{pidc} Kernel PCA (RBF) on top PCs')
        ax.grid(True, ls='--', lw=0.5, alpha=0.5)
        plt.show()

    # Return everything needed for downstream inspection and decisions:
    # aligned coordinates (for sanity checks), spectral quantities (eigvals/eigvecs),
    # model embeddings (scores/kpca_scores), and hinge profiles (wrmsd_modes).
    return {
        'aligned_coords': aligned,  # aligned coordinates (for sanity checks)
        'eigvals': eigvals,  # spectral quantities (eigvals/eigvecs)
        'eigvecs': eigvecs,  # spectral quantities (eigvals/eigvecs)
        'scores': scores,  # shape: (modes, M)
        'wrmsd_modes': wr,  # dict: mode -> per-residue amplitudes (Å)
        'kpca_scores': kpca_scores  # shape: (2, M) or None
    }


if __name__ == '__main__':

    rp_pidc_parsed_cif_dir = _rp_parsed_cifs('multimod_2713_hetallchains_hom1chain')
    rp_pidc_ssvs = sorted(glob.glob(os.path.join(rp_pidc_parsed_cif_dir, '*.ssv')))
    for rp_pidc_ssv in  rp_pidc_ssvs:
        pidc = os.path.basename(rp_pidc_ssv).removesuffix('.ssv')
        print(pidc)
        _pidc_pdf = pd.read_csv(rp_pidc_ssv, sep=' ')

        res = essential_dynamics_pca(
            pidc=pidc,
            pidc_pdf=_pidc_pdf,
            use_correlation=False,  # try both True/False; paper suggests comparing
            top_modes_for_rmsd=(0, 1, 2),
            try_kmeans_k=2,  # try 2 clusters; adjust as needed
            do_kpca=True  # turn on kernel PCA if clusters are “almost separated”
        )
