"""
Coarse-grained normal mode analysis (NMA) using the Gaussian network model (GNM).
- Builds Kirchhoff (Laplacian) matrix from alpha-carbon coordinates.
- Computes normal-mode decomposition of Kirchhoff matrix.
- Computes mean-square fluctuations (MSFs) and cross-correlations from Kirchhoff (Laplacian) matrix^{-1}.
- Plots profiles and matrices (helpers).

References (GNM / ENM theory):
- Bahar I., Atilgan AR., Erman B. Fold Des. 1997, 2, 173–181.  (original GNM)
- Atilgan AR. et al., Biophys J. 2001, 80, 505–515.           (ANM / ENM)
- Bahar I., Lezon TR., Yang LW., Eyal E. Annu Rev Biophys. 2010, 39, 23–42. “Gaussian network model” article (for formulas for MSF & correlations).

-------------------------------------------------------
Approximate mapping of Python functions to literature:
-------------------------------------------------------

- distance_matrix(), build_kirchhoff_matrix():
    Bahar et al., Fold Des 1997.
    Definition of the GNM Kirchhoff (Laplacian) matrix (Gamma) from alpha-carbon-based contact network.

- gnm_modes():
    Bahar et al., Fold Des 1997; Bahar et al., Annu Rev Biophys 2010.
    Eigen-decomposition of Kirchhoff matrix (Gamma), removal of translational zero mode, interpretation of
    low-frequency eigenmodes as collective motions.

- gnm_fluctuations():
    Bahar et al., Fold Des 1997; Bahar et al., Annu Rev Biophys 2010.
    Mean-square fluctuations: Gamma^{-1} = U Λ^{+} U-transposed,
    Cross-correlations: <ΔR_i^2> ∝ (Gamma^{-1})_{ii}
    And their normalisation: <ΔR_i·ΔR_j> ∝ (Gamma^{-1})_{ij},

- msf_to_bfactors():
    Bahar et al., Fold Des 1997.
    Relation B_i = (8π^2 / 3) <ΔR_i^2>.

- run_gnm():
    Bahar et al., Annu Rev Biophys 2010.
    Overall workflow from structure to Gamma to modes to fluctuations, correlations, and B-factors.

- Atilgan et al., Biophys J 2001 (ANM/ENM):
    Not directly implemented here; would correspond to replacing scalar Gamma with a 3N×3N anisotropic Hessian
    and performing ANM rather than isotropic GNM.
"""

import os, glob
from typing import Optional, Tuple
from dataclasses import dataclass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _rp_nmr_dir() -> str:
    return os.path.join('..', 'data', 'NMR')

def _rp_parsed_cifs(subdir: str) -> str:
    return os.path.join(_rp_nmr_dir(), 'parsed_cifs', subdir)


# --------------------------------------------------------------------
# Core linear algebra
# --------------------------------------------------------------------

def _distance_matrix(coords: np.ndarray) -> np.ndarray:
    """
    Compute full pairwise Euclidean distance matrix for N nodes.

    :param coords: (N,3). Cartesian coordinates (alpha-carbons only).
    :return: (N,N). Distance matrix, D[i, j] = ||r_i - r_j||_2.
    """
    coords = np.asarray(coords, dtype=float)  # already done in `run_gnm()`, i.e. defensive coding.
    # Using None in this way is just another way (than e.g. using `.reshape()`) for reshaping coords (N,3) to diff (N,N,3).
    coords_col_vec = coords[:, None, :]  # (N,1,3) <- (N,3)
    coords_row_vec = coords[None, :, :]  # (1,N,3) <- (N,3)
    diff = coords_col_vec - coords_row_vec  # (N,N,3) <- (N,1,3) - (1,N,3)
    dist_matrx = np.linalg.norm(diff, axis=-1)  # (N,N) <- (N,N,3)
    np.fill_diagonal(dist_matrx, 0.0)  # should already have 0.0 along diagonal, but just reinforces this.
    return dist_matrx

def _build_kirchhoff_matrix(coords: np.ndarray, cutoff: float = 7.5, gamma: float = 1.0) -> np.ndarray:
    """
    (Note: Don't confuse the spring constant `gamma` (lowercase g) with the Laplacian matrix `Gamma` (uppercase G).)

    Build the Gaussian network model Kirchhoff (Laplacian) matrix:

        Gamma_ij = -gamma if i≠j and distance(i,j) < cutoff 0.0
        otherwise Gamma_ii = -∑_{j≠i} Gamma_ij   (degree on diagonal)

    :param coords: (N,3). Cartesian node positions (alpha-carbons only).
    :param cutoff: Distance cutoff (Å) for defining contacts / springs. Standard GNM values are ~7–8 Å for
    residue-level networks.
    :param gamma: Spring constant (units cancel in *relative* fluctuations, but kept for completeness).
    :return: (N,N). Kirchhoff / graph Laplacian matrix (Gamma).
    """
    # Compute the adjaceny matrix:
    D = _distance_matrix(coords)
    N = D.shape[0]
    contact = (D < cutoff) & (D > 0.0)  # Contact / adjacency mask (exclude self-contacts)

    # Compute the weighted Laplacian by: Gamma = gamma(degree matrix - adjacency matrix)
    Gamma = np.zeros((N, N), dtype=float)
    Gamma[contact] = -gamma  # Off-diagonals: -γ where there is a contact
    np.fill_diagonal(Gamma, -Gamma.sum(axis=1))  # Diagonals: degree = -sum(off-diagonals)
    return Gamma

@dataclass
class GNMResults:
    """Container for GNM outputs."""
    coords: np.ndarray               # (N,3) alpha-carbon coordinates
    Gamma: np.ndarray                # (N,N) Kirchhoff matrix (not same as spring constant lowercase `gamma`)
    eigvals: np.ndarray              # (N-1,) non-zero eigenvalues (ascending)
    eigvecs: np.ndarray              # (N,N-1) eigenvectors (columns) for non-zero modes
    msf: np.ndarray                  # (N,) mean-square fluctuations
    corr: np.ndarray                 # (N,N) cross-correlations (normalised)
    B_factors: np.ndarray            # (N,) predicted B-factors (8π^2/3 * msf)


# --------------------------------------------------------------------
# Normal modes & fluctuations
# --------------------------------------------------------------------

def _gnm_modes(Gamma: np.ndarray, n_modes: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    (Note: Spring constant `gamma` is with lowercase 'g'. Laplician matrix `Gamma` is with uppercase 'G'.)

    Diagonalise the Kirchhoff matrix and return non-zero modes.
    (A zero mode is an eigenvector whose eigenvalue is exactly zero. It does not represent internal motions (like
    non-zero modes) but rather rigid-body motions, i.e. rotations/translations of the whole protein as one, so they're
    global symmetries.
    The next function called, `_gnm_fluctuations()`, computes the pseudo-inverse, by eigendecomposition, which would
    implicitly remove zero modes. However, the zero modes are still explicitly removed in this function, for
    convenience and to avoid any division by zero.

    :param Gamma: (N,N). Kirchhoff matrix.
    :param n_modes: Number of *non-zero* modes to keep (lowest-frequency ones). If None, return all N-1 non-zero modes.
    :return:
            eigvals_nz : (K,). Non-zero eigenvalues (ascending; small -> soft, collective modes).
            eigvecs_nz : (N,K). Corresponding eigenvectors (columns). Mode k is eigvecs_nz[:, k].
    """
    Gamma = np.asarray(Gamma, dtype=float)
    eigvals, eigvecs = np.linalg.eigh(Gamma)  # Symmetric positive semi-definite
    order = np.argsort(eigvals)  # Sort ascending (for GNM, smallest non-zero eigenvalues = softest modes)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    # The first eigenvalue should be ~0 (global translation)
    # Drop zero (or numerically ~zero) modes.
    tol = 1e-6
    non_zero = eigvals > tol
    eigvals_nz = eigvals[non_zero]
    eigvecs_nz = eigvecs[:, non_zero]

    if n_modes is not None:
        n = min(n_modes, eigvals_nz.shape[0])
        eigvals_nz = eigvals_nz[:n]
        eigvecs_nz = eigvecs_nz[:, :n]

    return eigvals_nz, eigvecs_nz


def _gnm_fluctuations(eigvals: np.ndarray, eigvecs: np.ndarray, kBT: float = 1.0, gamma: float = 1.0,
                      normalise_corr: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    (Note: Spring constant `gamma` is with lowercase 'g'. Laplician matrix `Gamma` is with uppercase 'G'.)

    Compute mean square fluctuations (MSFs) and cross-correlations from GNM eigenpairs.
    Convert the non-zero normal modes of a GNM Laplacian into predicted equilibrium residue fluctuations and
    correlations by constructing the pseudo-inverse of the Kirchhoff matrix and interpreting it as a covariance
    operator under a harmonic, locally stable model of the protein.

    Using:
        Gamma^{-1} = U Λ^{+} U^T  (pseudo-inverse, zero mode excluded)
            where `Λ^{+}` is `diag(1 / λ_k)` and `λ_k` is the k-th eigenvalue,
            and `U^T` is `U` transposed.

        <ΔR_i · ΔR_j> = (3 k_B T / gamma) (Gamma^{-1})_{ij}
        <ΔR_i^2>      = (3 k_B T / gamma) (Gamma^{-1})_{ii}
            where `R_i` is the Cartesian position vector of residue `i`,
            `k_B T` is the thermal energy, `gamma` is the spring constant,
            `Gamma^{-1}` is the pseudo-inverse of the Laplacian matrix,
            `ΔR_i` is the displacement of residue `i` from its equilibrium position,
            `<ΔR_i · ΔR_j>` is the average correlated motion between residues `i` and `j`
            `<ΔR_i^2>` is the average squared fluctuation amplitude of residue `i`, averaged over all thermally
            accessible fluctuations. This is called the mean-squared fluctions (MSF).


    :param eigvals: Non-zero eigenvalues (ascending), shape (K,).
    :param eigvecs: Corresponding eigenvectors (columns), shape (N,K).
    :param kBT: Thermal energy (k_B T). Often just set to 1 for relative values.
    :param gamma: Spring constant used in building Gamma.
    :param normalise_corr: If True, return correlation matrix normalised to [-1, 1]:
                            C_ij = <ΔR_i · ΔR_j> / sqrt( <ΔR_i^2> <ΔR_j^2> ).
    :return:
        msf: The predicted mean-square fluctuations (MSF), <ΔR_i^2>, (shape (N,)) of each residue i around the assumed
        equilibrium structure.
        cov: The full covariance matrix, <ΔR_i · ΔR_j> (shape (N,N)) of residue displacements predicted by the harmonic
        GNM model.
        corr: Information about correlated residue motions: either the physically scaled covariance matrix (when
        `normalise_corr=False`, where amplitude and directionality are coupled) and hence duplicating the other
        returned value `cov`, or the normalised cross-correlation matrix (when `normalise_corr=True`), which isolates
        relative directionality in the range [−1,1].
    """
    eigvals = np.asarray(eigvals, dtype=float)
    eigvecs = np.asarray(eigvecs, dtype=float)

    # Pseudo-inverse of Gamma (Γ) via eigen-decomposition Γ^{-1} = U Λ^{+} U^T with Λ^{+}_kk = 1/λ_k
    # is done in the following 2 lines:
    inv_lambda = 1.0 / eigvals  # 1st: compute Λ^{+}, the diagonal of inverted eigenvalues for non-zero modes, (shape (K,)).
    # cov0 = U diag(Λ^{+}) U^T  # cov0 covariance kernel, before physical scaling.
    cov0 = eigvecs @ (inv_lambda[:, None] * eigvecs.T)  # 2nd: scale eigenvectors by Λ^{+} and reconstruct Gamma^{-1}.

    pref = 3.0 * kBT / gamma  # Physical prefactor; for relative fluctuations, you can set kBT/gamma = 1.
    cov = pref * cov0  # Physically scaled covariance.
    msf = np.diag(cov).copy()  # Mean-square fluctuations = diagonal elements

    corr = cov
    if normalise_corr:
        # Normalise to get correlation coefficients in [-1,1]
        std = np.sqrt(msf)
        std[std == 0] = 1.0  # Avoid division by zero
        outer_std = np.outer(std, std)
        corr = corr / outer_std
        np.clip(corr, -1.0, 1.0, out=corr)  # Numerical clipping
    return msf, cov, corr


def _msf_to_bfactors(msf: np.ndarray) -> np.ndarray:
    """
    Convert mean-square fluctuations to equivalent isotropic B-factors via:
        B_i = (8π^2 / 3) <ΔR_i^2>

    :param msf: (N,). Mean-square fluctuations (columns).
    :return: (N,). Predicted equivalent isotropic B-factors.
    """
    return (8.0 * np.pi**2 / 3.0) * msf


# --------------------------------------------------------------------
# Plotting helpers
# --------------------------------------------------------------------

def _plot_msf_and_bfactors(msf: np.ndarray, title: str = "") -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    x = np.arange(1, len(msf) + 1)
    factor = (8.0 * np.pi**2) / 3.0

    fig, ax1 = plt.subplots()
    ax1.plot(x, msf)
    ax1.set_xlabel("Residue index")
    ax1.set_ylabel(r'<$\Delta R_i^2$> (arb. units)')
    ax1.set_title(title)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.set_ylabel('B (Å$^2$, arbitrary scale)')

    # Match B-factor scale to MSF scale
    y1_min, y1_max = ax1.get_ylim()
    ax2.set_ylim(y1_min * factor, y1_max * factor)

    plt.show()


def _plot_sequence_signal(signal: np.ndarray, title: str, ylabel: str = 'value'):
    """Simple per-residue line plot."""
    signal = np.asarray(signal)
    fig, ax = plt.subplots(figsize=(10, 2.8))
    ax.plot(np.arange(1, signal.shape[0] + 1), signal, lw=1.5)
    ax.set_xlabel('Residue index')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.margins(x=0)
    ax.grid(ls='--', lw=0.5, alpha=0.5)
    plt.show()
    return fig, ax


def _plot_matrix(M: np.ndarray, title: str, cmap: str = 'RdBu_r', vmin=None, vmax=None):
    """Heatmap of an N×N matrix (e.g. correlation matrix)."""
    M = np.asarray(M)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(M, origin='upper', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Residue index')
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.show()
    return fig, ax


# --------------------------------------------------------------------
# Main function
# --------------------------------------------------------------------

def run_gnm(coords: np.ndarray, cutoff: float = 7.5, gamma: float = 1.0, n_modes: Optional[int] = None,
            kBT: float = 1.0, do_plots: bool = True, pidc: Optional[str] = None) -> GNMResults:
    """
    (Note: Spring constant `gamma` is with lowercase 'g'. Laplician matrix `Gamma` is with uppercase 'G'.)

    High-level driver: from alpha-carbon coordinates → GNM modes & fluctuation profiles.

    :param coords: (NxM,3). Alpha-carbon coordinates for a single structure (one model, one chain).
    :param cutoff: Distance cutoff (Å) for defining springs in the network.
    :param gamma: Spring constant in the Kirchhoff matrix.
    :param n_modes: Number of non-zero modes to retain (lowest-frequency). If None, use all N-1 non-zero modes.
    :param kBT: Thermal energy; for relative profiles, 1.0 is fine.
    :param do_plots: True to plot MSF profile and correlation matrix.
    :param pidc: Optional identifier for titles (e.g. '1A03_A').
    :return: GNMResults. Dataclass with all relevant outputs.
        coords:     (N,3) alpha-carbon coordinates (given as input argument to this function).
        Gamma:      (N,N) Kirchhoff matrix
        eigvals:    (N-1,) non-zero eigenvalues
        eigvecs:    (N,N-1) eigenvectors (columns) for non-zero modes
        msf:        (N,) mean-square fluctuations
        corr:       (N,N) cross-correlations (normalised)
        B_factors:  (N,) predicted B-factors (8π^2/3 * msf)
    """
    coords = np.asarray(coords, dtype=float)
    Gamma = _build_kirchhoff_matrix(coords, cutoff=cutoff, gamma=gamma)
    eigvals, eigvecs = _gnm_modes(Gamma, n_modes=n_modes)  # Diagonalise and get non-zero modes
    msf, cov, corr = _gnm_fluctuations(eigvals, eigvecs, kBT=kBT, gamma=gamma)  # fluctuations and correlations
    bfacs = _msf_to_bfactors(msf)  # Equivalent B-factors

    if do_plots:
        # _plot_sequence_signal(msf, title=f'GNM mean-square fluctuations for {pidc}', ylabel=r'<$\Delta R_i^2$> (arb. units)')
        # _plot_sequence_signal(bfacs, title=f'GNM predicted B-factors for {pidc}', ylabel='B (Å², arbitrary scale)')
        _plot_msf_and_bfactors(msf, title=f'GNM predicted MSFs & B-factors for {pidc}')
        _plot_matrix(corr, title=f'GNM cross-correlations for {pidc}', cmap="RdBu_r", vmin=-1.0, vmax=1.0)

    return GNMResults(coords=coords, Gamma=Gamma, eigvals=eigvals, eigvecs=eigvecs, msf=msf, corr=corr, B_factors=bfacs)


if __name__ == '__main__':

    rp_pidc_parsed_cif_dir = _rp_parsed_cifs('multimod_2713_hetallchains_hom1chain')
    rp_pidc_ssvs = sorted(glob.glob(os.path.join(rp_pidc_parsed_cif_dir, '*.ssv')))
    cols = ['A_Cartn_x', 'A_Cartn_y', 'A_Cartn_z']
    model_num = 1

    for rp_pidc_ssv in rp_pidc_ssvs[2:3]:
        pidc = os.path.basename(rp_pidc_ssv).removesuffix('.ssv')
        print(pidc)
        pidc_pdf = pd.read_csv(rp_pidc_ssv, sep=' ')
        pidc_pdf = pidc_pdf[pidc_pdf['A_pdbx_PDB_model_num'] == model_num]
        pidc_coords = pidc_pdf[cols].to_numpy(dtype=float)
        results = run_gnm(coords=pidc_coords, pidc=pidc)
