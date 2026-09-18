import numpy as np
import matplotlib.pyplot as plt


def compute_kernel_matrix_vectorized(r_h_vec, r_p_vec):
    """Calcule la matrice noyau K_ij de façon vectorisée."""
    r_h = np.asarray(r_h_vec)[:, np.newaxis]
    r_p = np.asarray(r_p_vec)[np.newaxis, :]
    return np.maximum(0, 1.0 - (r_h / r_p)) ** 2


def deconvolve_richardson_lucy(Ksec_exp, K_phys, iterations=500):
    """
    Algorithme itératif de Richardson-Lucy.
    La régularisation est contrôlée par le nombre d'itérations.
    """
    n_p = K_phys.shape[1]

    # 1. Initialisation : On part d'une distribution totalement plate (uniforme)
    psd = np.ones(n_p) / n_p

    # 2. Pré-calculs pour accélérer la boucle
    K_T = K_phys.T
    eps = 1e-12  # Pour éviter les divisions par zéro

    # Facteur de normalisation de l'opérateur K
    K_norm = np.sum(K_phys, axis=0)
    K_norm[K_norm == 0] = eps

    # Historique de l'erreur pour visualiser la convergence
    errors = []

    # 3. Boucle itérative
    for _ in range(iterations):
        # Modèle direct : Ksec théorique généré par la PSD actuelle
        Ksec_est = K_phys @ psd

        # Calcul de l'erreur (Norme L2)
        errors.append(np.linalg.norm(Ksec_est - Ksec_exp))

        # Ratio entre la cible expérimentale et l'estimation actuelle
        ratio = Ksec_exp / (Ksec_est + eps)

        # Mise à jour multiplicative de la PSD
        psd = psd * (K_T @ ratio) / K_norm

    return psd, errors


# ==============================================================================
# EXECUTION DU SCRIPT
# ==============================================================================
if __name__ == "__main__":
    np.random.seed(42)

    # 1. Préparation des données synthétiques
    r_h_exp = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0])
    r_p_vec = np.linspace(0.6, 50.0, 500)
    dr_p = r_p_vec[1] - r_p_vec[0]

    # Vraie PSD bimodale (purement pour vérification visuelle)
    true_psd = 0.6 * np.exp(-(((r_p_vec - 5.0) / 2.0) ** 2)) + 0.4 * np.exp(
        -(((r_p_vec - 25.0) / 6.0) ** 2)
    )
    true_psd /= np.trapz(true_psd, r_p_vec)

    K_phys = compute_kernel_matrix_vectorized(r_h_exp, r_p_vec)
    Ksec_pure = K_phys @ (true_psd * dr_p)

    # Ajout d'un bruit réaliste (2%)
    noise = np.random.normal(0, 0.02, size=len(r_h_exp))
    Ksec_exp = np.clip(Ksec_pure + noise, 0.0, 1.0)

    # 2. Déconvolution par Richardson-Lucy
    # Tu peux modifier ce chiffre pour voir l'effet (ex: 50, 500, 5000)
    n_iterations = 800
    w_rec, erreurs = deconvolve_richardson_lucy(
        Ksec_exp, K_phys, iterations=n_iterations
    )

    # Recalcul de la courbe ajustée pour l'affichage
    Ksec_fit = K_phys @ w_rec
    psd_density = w_rec / (np.trapz(w_rec, r_p_vec) + 1e-12)

    # 3. Affichage des résultats
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 5))

    ax1.plot(r_h_exp, Ksec_exp, "ro", label="Ksec exp (bruité)")
    ax1.plot(r_h_exp, Ksec_fit, "b--", label=f"Ajustement RL")
    ax1.set_title("Exclusion ISEC")
    ax1.legend()

    ax2.plot(r_p_vec, true_psd, "k-", label="Vraie PSD")
    ax2.plot(
        r_p_vec, psd_density, "g-", linewidth=2, label=f"PSD RL ({n_iterations} iters)"
    )
    ax2.set_title("Distribution de Taille des Pores")
    ax2.legend()

    plt.tight_layout()
    plt.show()
