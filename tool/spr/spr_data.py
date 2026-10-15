from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Source des fenetres et des libelles :
# https://raw.githubusercontent.com/jejessika/webappspr/main/190520_Hasil%20Algo.csv
PHASES = [
    {"echantillon": "RBD", "t_debut_s": 1280., "t_fin_s": 2459., "concentration": "14 ng/mL"},
    {"echantillon": "IBV", "t_debut_s": 3361., "t_fin_s": 4508., "concentration": "3 ng/mL"},
    {"echantillon": "Non specifique", "t_debut_s": 5340., "t_fin_s": 6520., "concentration": "15 ng/mL"},
    {"echantillon": "Inconnu", "t_debut_s": 7154., "t_fin_s": 8530., "concentration": "0 ng/mL (Negatif)"},
]

BASELINE_OFFSET_MIN = 120  # s avant la borne gauche de la fenetre
BASELINE_OFFSET_MAX = 300  # s avant la borne gauche de la fenetre
PLATEAU_DUREE = 120       # derniere tranche de la fenetre

def charger_csv(chemin):
    # index_col=False est necessaire : le fichier comporte une virgule finale
    # supplementaire sur chaque ligne; usecols seul decalera les valeurs.
    donnees = pd.read_csv(chemin, index_col=False)
    x = pd.to_numeric(donnees.iloc[:, 0], errors="coerce")
    y = pd.to_numeric(donnees.iloc[:, 1], errors="coerce")
    tableau = pd.DataFrame({"temps_s": x, "reponse_RU": y})
    tableau = tableau.replace([np.inf, -np.inf], np.nan).dropna()
    tableau = tableau.sort_values("temps_s").drop_duplicates("temps_s")
    return tableau.reset_index(drop=True)


def analyser(df, phases):
    t = df["temps_s"].to_numpy()
    y = df["reponse_RU"].to_numpy()
    resume = []
    profils = []
    for phase in phases:
        nom = phase["echantillon"]
        a, b = phase["t_debut_s"], phase["t_fin_s"]
        avant = (t >= a - BASELINE_OFFSET_MAX) & (t < a - BASELINE_OFFSET_MIN)
        fin = (t >= b - PLATEAU_DUREE) & (t <= b)
        pendant = (t >= a) & (t <= b)
        apres = (t > b + 120) & (t <= b + 300)
        if not (avant.any() and fin.any() and pendant.any()):
            print(f"Fenetre incomplete ignoree: {nom}")
            continue

        niveau_avant = float(np.median(y[avant]))
        niveau_fin = float(np.median(y[fin]))
        maximum = float(np.max(y[pendant]))
        niveau_apres = float(np.median(y[apres])) if apres.any() else np.nan
        resume.append({
            "echantillon": nom,
            "concentration_annotee": phase["concentration"],
            "debut_s": a,
            "fin_s": b,
            "RU_avant_mediane": niveau_avant,
            "RU_fin_mediane": niveau_fin,
            "delta_RU_fin_moins_avant": niveau_fin - niveau_avant,
            "maximum_RU_pendant": maximum,
            "RU_apres_mediane": niveau_apres,
            "retour_apres_moins_avant_RU": niveau_apres - niveau_avant,
        })
        profils.append((phase, t[pendant] - a, y[pendant] - niveau_avant))
    return pd.DataFrame(resume), profils


def tracer(df, phases, profils, fichier_png):
    t = df["temps_s"].to_numpy()
    y = df["reponse_RU"].to_numpy()
    couleurs = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9))

    ax1.plot(t, y, linewidth=1.5, color="black", label="Signal traite (RU)")
    for i, phase in enumerate(phases):
        a, b = phase["t_debut_s"], phase["t_fin_s"]
        couleur = couleurs[i % len(couleurs)]
        label = f'{phase["echantillon"]} - {phase["concentration"]}'
        ax1.axvspan(a, b, alpha=0.15, color=couleur, label=label)
        ax1.axvline(a, linewidth=0.8, linestyle=":", color=couleur)
        ax1.axvline(b, linewidth=0.8, linestyle=":", color=couleur)

    ax1.set(title="Sensorigramme SPR annote - fichier webappspr", xlabel="Temps (s)", ylabel="Reponse (RU)")
    ax1.grid(alpha=0.2)
    ax1.legend(fontsize=9, loc="upper left", ncol=2)

    for i, (phase, temps_relatif, y_relatif) in enumerate(profils):
        ax2.plot(temps_relatif, y_relatif, linewidth=1.6,
                 color=couleurs[i % len(couleurs)],
                 label=f'{phase["echantillon"]} - {phase["concentration"]}')
    ax2.axhline(0, linestyle="--", linewidth=0.8, color="black")
    ax2.set(title="Reponses dans les fenetres annotees (baseline locale soustraite)",
            xlabel="Temps depuis le debut de la fenetre (s)",
            ylabel="Reponse relative a la baseline locale (RU)")
    ax2.grid(alpha=0.2)
    ax2.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(fichier_png, dpi=180)
    plt.show()


chemin = Path("190520_.csv")
dossier_sortie = chemin.parent
df = charger_csv(chemin)
resume, profils = analyser(df, PHASES)
df.to_csv(dossier_sortie / "spr_temps_ru_nettoye.csv", index=False)
resume.to_csv(dossier_sortie / "spr_resume_echantillons.csv", index=False)
tracer(df, PHASES, profils, dossier_sortie / "spr_sensorigramme_annote.png")
print(f"Points conserves : {len(df)}")
print(resume.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
print("\nLes variations sont descriptives; elles ne sont pas des mesures de K_D, k_a ou k_d.")


# %%

