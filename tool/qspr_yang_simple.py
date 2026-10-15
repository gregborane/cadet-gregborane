#!/usr/bin/env python3
"""QSPR–SVR de ligands multimodes, version scientifique fonctionnelle.

Installation :
    python -m pip install numpy pandas scipy scikit-learn rdkit matplotlib
    python -m pip install mordredcommunity  # seulement si BACKEND = 'mordred'
Exécution : python qspr_yang_simple.py

Changer les paramètres ci-dessous. CSV : smiles,target et éventuellement ph,
ligand_id,protonated_smiles. Une réponse/protéine/protocole par modèle.
Pour les ligands réels, inclure le fragment de linker comme dans le papier,
et fournir de préférence les états protonés préparés à partir de pKa adaptés.

Référence : Yang, Breneman, Cramer (2007), J. Chromatogr. A 1175, 96–105,
DOI 10.1016/j.chroma.2007.10.037, sections 2.1–2.3, pp. 97–99.
Le papier utilise 132 descripteurs MOE, MMFF94X et ACD/pKa : les substituts
ouverts ci-dessous ne reproduisent pas exactement cette banque. Le rayon
MOE discuté p. 104 est topologique; le rayon de giration ajouté ici est 3D.
Le papier nomme 'bootstrapping' des combinaisons train/validation : notre
mode 'subsample' en est une approximation aléatoire, non exhaustive.
Le mode 'bootstrap' réalise un tirage classique avec remise.
Pas de classes personnalisées, de try/except ou de contrôles d'erreurs.
Les erreurs de données ou de calcul remontent directement des bibliothèques.
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Descriptors3D
from scipy import sparse
from scipy.optimize import linprog
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ------------------------- Paramètres scientifiques -------------------------
CSV = None                         # ex. 'mes_ligands.csv'; None = jeu fictif
OUTPUT = Path('resultats_qspr_simple')
BACKEND = 'rdkit'                   # 'rdkit' ou 'mordred'
PH = 6.8
N_RESAMPLES = 50
N_FOLDS = 5
TRAIN_FRACTION = 0.8
EXTERNAL_FRACTION = 0.2
RESAMPLING = 'subsample'            # 'subsample' ou 'bootstrap'
SEED = 42
L1_C = 1.0                         # objectif : ||w||1 + C * somme des pertes
L1_EPSILON = 0.1                   # epsilon en unités de y standardisé
SVR_C = 10.0
SVR_EPSILON = 0.1                   # epsilon en unités de y standardisé
SVR_GAMMA = 'scale'
COEF_TOL = 1e-7                     # seuil numérique des coefficients non nuls
UNIT = 'mS/cm'


def toy_data():
    """15 ligands simples; cibles entièrement fictives, aucune donnée du papier."""
    smiles = ['CC(=O)O', 'CCC(=O)O', 'CCCC(=O)O', 'CCCCC(=O)O',
              'O=C(O)c1ccccc1', 'O=C(O)Cc1ccccc1', 'O=C(O)CCc1ccccc1',
              'OCC(=O)O', 'O=C(O)CCO', 'NCC(=O)O', 'CC(N)C(=O)O',
              'NCCC(=O)O', 'O=C(O)Cc1ccc(O)cc1', 'O=C(O)c1ccncc1',
              'O=C(O)C1CCCCC1']
    return pd.DataFrame({'smiles': smiles, 'ph': PH,
                         'target': [11, 15, 19, 24, 29, 33, 38, 9, 12,
                                    17, 20, 22, 27, 21, 30]})


def protonate(smiles, ph):
    """État dominant approximatif; pKa génériques, pas prédiction de pKa.

    Carboxyle 4.5, phénol 10, amine aliphatique 9.5, pyridine 5.2.
    Pas de couplage de sites ni de spéciation complète. Utiliser la colonne
    protonated_smiles pour remplacer ces règles par des états validés.
    """
    mol = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
    rw = Chem.RWMol(mol)
    acids = [('[CX3](=O)[O;H1,-1]', 2, 4.5), ('[c][O;H1,-1]', 1, 10.)]
    for pattern, position, pka in acids:
        for match in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
            atom = rw.GetAtomWithIdx(match[position])
            atom.SetFormalCharge(-int(ph >= pka))
            atom.SetNumExplicitHs(int(ph < pka))
            atom.SetNoImplicit(True)
    bases = [('[N;X3,X4;+0,+1;!$(N-C=O);!$(N-S(=O)=O);!$(N-a)]', 9.5),
             ('[nH0+0,nH+1;r6]', 5.2)]
    for pattern, pka in bases:
        for match in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
            old = mol.GetAtomWithIdx(match[0])
            atom = rw.GetAtomWithIdx(match[0])
            neutral_h = old.GetTotalNumHs() - int(old.GetFormalCharge() == 1)
            atom.SetFormalCharge(int(ph < pka))
            atom.SetNumExplicitHs(neutral_h + int(ph < pka))
            atom.SetNoImplicit(True)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def molecular_descriptors(data, backend='rdkit'):
    """Banque 2D/3D, Gasteiger, PEOE_VSA, SlogP et indices topologiques.

    Un conformère ETKDG puis MMFF94 (UFF si paramètres MMFF indisponibles).
    Les objets Missing/Error de Mordred deviennent NaN avant nettoyage.
    """
    calculator = None
    if backend == 'mordred':
        from mordred import Calculator, descriptors
        calculator = Calculator(descriptors, ignore_3D=False)
    rows, audit = [], []
    for _, row in data.iterrows():
        ph = row.get('ph', PH)
        provided = row.get('protonated_smiles', None)
        mol = (Chem.MolFromSmiles(provided) if pd.notna(provided)
               else protonate(row.smiles, ph))
        values = {'RD_' + name: fun(mol) for name, fun in Descriptors.descList}
        mh = Chem.AddHs(mol)
        AllChem.ComputeGasteigerCharges(mh)
        q = np.array([float(a.GetProp('_GasteigerCharge')) for a in mh.GetAtoms()])
        values.update(Gasteiger_min=q.min(), Gasteiger_max=q.max(),
                      Gasteiger_std=q.std(), Gasteiger_abs_sum=np.abs(q).sum())
        params = AllChem.ETKDGv3()
        params.randomSeed = SEED
        AllChem.EmbedMolecule(mh, params)
        if AllChem.MMFFHasAllMoleculeParams(mh):
            AllChem.MMFFOptimizeMolecule(mh, maxIters=1000)
        else:
            AllChem.UFFOptimizeMolecule(mh, maxIters=1000)
        values.update({'3D_' + name: fun(mh) for name, fun in Descriptors3D.descList})
        if calculator is not None:
            result = calculator(mh)
            values.update({'M_' + str(d): v for d, v in zip(calculator.descriptors, result)})
        values['condition_pH'] = ph
        rows.append(values)
        audit.append({'smiles': row.smiles, 'ph': ph,
                      'used_smiles': Chem.MolToSmiles(mol),
                      'formal_charge': Chem.GetFormalCharge(mol)})
    X = pd.DataFrame(rows).apply(pd.to_numeric, errors='coerce')
    return X.replace([np.inf, -np.inf], np.nan), pd.DataFrame(audit)


def clean_training(X):
    """Élimination des colonnes NaN et constantes sur entraînement seulement."""
    Z = X.dropna(axis=1)
    Z = Z.loc[:, Z.var(ddof=0) > 1e-12]
    return Z, Z.mean()


def resample_indices(n, rng, mode):
    """Sous-échantillonnage sans remise (papier) ou bootstrap avec remise."""
    if mode == 'bootstrap':
        return rng.integers(0, n, n)
    return rng.choice(n, size=int(TRAIN_FRACTION*n), replace=False)


def linear_l1_svr(X, y, C=L1_C, epsilon=L1_EPSILON):
    """SVR linéaire ℓ1 résolu par LP : ||w||1 + C Σ(ξ+ξ*).

    w = w+ - w-, intercept libre; contrairement à LinearSVR, les
    coefficients sont pénalisés ℓ1. La perte est epsilon-insensible.
    """
    n, p = X.shape
    sx = sparse.csr_matrix(X)
    one = sparse.csr_matrix(np.ones((n, 1)))
    eye = sparse.eye(n, format='csr')
    zero = sparse.csr_matrix((n, n))
    A = sparse.vstack([sparse.hstack([-sx, sx, -one, -eye, zero]),
                       sparse.hstack([sx, -sx, one, zero, -eye])], format='csr')
    objective = np.r_[np.ones(2*p), 0., np.full(2*n, C)]
    solution = linprog(objective, A_ub=A, b_ub=np.r_[epsilon-y, epsilon+y],
                       bounds=[(0, None)]*(2*p) + [(None, None)] + [(0, None)]*(2*n),
                       method='highs')
    return solution.x[:p] - solution.x[p:2*p]


def select_union(X, y, rng, n_resamples, mode):
    """UNION de tous les coefficients non nuls; conserve les poids par tirage.

    Les observations exclues de chaque tirage ne participent ni aux scalers
    ni à l'ajustement. Pas de sélection par fréquence ou de top-k.
    """
    weights = pd.DataFrame(0., index=range(n_resamples), columns=X.columns)
    for iteration in range(n_resamples):
        idx = resample_indices(len(y), rng, mode)
        Z, _ = clean_training(X.iloc[idx])
        xb = StandardScaler().fit_transform(Z)
        yb = StandardScaler().fit_transform(y[idx, None]).ravel()
        weights.loc[iteration, Z.columns] = linear_l1_svr(xb, yb)
    selected = weights.columns[(weights.abs() > COEF_TOL).any(axis=0)].tolist()
    return selected, weights


def fit_qspr(X, y, seed=SEED, n_resamples=N_RESAMPLES, mode=RESAMPLING):
    """Nettoyage, sélection UNION puis ensemble SVR RBF; scalers par modèle."""
    rng = np.random.default_rng(seed)
    Z, means = clean_training(X)
    selected, weights = select_union(Z, y, rng, n_resamples, mode)
    Z = Z[selected]
    models = []
    for _ in range(n_resamples):
        idx = resample_indices(len(y), rng, mode)
        sx = StandardScaler().fit(Z.iloc[idx])
        sy = StandardScaler().fit(y[idx, None])
        svr = SVR(kernel='rbf', C=SVR_C, epsilon=SVR_EPSILON, gamma=SVR_GAMMA)
        svr.fit(sx.transform(Z.iloc[idx]), sy.transform(y[idx, None]).ravel())
        models.append((sx, sy, svr))
    return {'selected': selected, 'means': means[selected],
            'weights': weights, 'models': models}


def predict_qspr(model, X):
    """Moyenne et dispersion entre modèles dans les unités originales.

    Pour les NaN de nouveaux ligands : moyenne d'entraînement.
    L'écart-type n'est pas un intervalle de prédiction calibré.
    """
    Z = X[model['selected']].fillna(model['means'])
    predictions = np.array([sy.inverse_transform(svr.predict(sx.transform(Z))[:, None]).ravel()
                            for sx, sy, svr in model['models']])
    return predictions.mean(axis=0), predictions.std(axis=0, ddof=1)


def cross_validation(X, y, groups, n_folds=N_FOLDS):
    """Tout le pipeline réappris à chaque fold, sans sélection globale préalable."""
    predicted, std = np.empty(len(y)), np.empty(len(y))
    for fold, (train, test) in enumerate(GroupKFold(n_folds).split(X, y, groups)):
        model = fit_qspr(X.iloc[train], y[train], seed=SEED+fold+1)
        predicted[test], std[test] = predict_qspr(model, X.iloc[test])
    return predicted, std


def scores(y, predicted):
    """R² global sur prédictions hors fold; RMSE dans les unités de la cible."""
    return {'R2': r2_score(y, predicted),
            'RMSE': np.sqrt(mean_squared_error(y, predicted))}


def plot_results(y, predicted, std, importance, metrics, output):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.errorbar(y, predicted, yerr=std, fmt='o', capsize=3)
    bounds = [min(y.min(), (predicted-std).min()), max(y.max(), (predicted+std).max())]
    ax.plot(bounds, bounds, '--', color='gray')
    ax.set(xlabel=f'Expérimental ({UNIT})', ylabel=f'Prédit CV ({UNIT})',
           title=f"R² CV = {metrics['R2']:.3f}; RMSE = {metrics['RMSE']:.2f}\nBarres : ±1 écart-type des modèles")
    fig.tight_layout()
    fig.savefig(output/'correlation_cv.png', dpi=180)
    plt.close(fig)
    top = importance.head(20).sort_values('relative_weight')
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(top.index, 100*top.relative_weight)
    ax.set(xlabel='Poids absolu moyen relatif (%)',
           title='Importance du sélecteur linéaire ℓ1 (20 premiers)')
    fig.tight_layout()
    fig.savefig(output/'importance_l1.png', dpi=180)
    plt.close(fig)


def main():
    data = pd.read_csv(CSV) if CSV is not None else toy_data()
    X, audit = molecular_descriptors(data, BACKEND)
    y = data.target.to_numpy(dtype=float)
    groups = (data.ligand_id.to_numpy() if 'ligand_id' in data.columns
              else np.array([Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in data.smiles]))
    # Test externe 'true unknowns' séparé avant toute sélection de variables.
    split = GroupShuffleSplit(n_splits=1, test_size=EXTERNAL_FRACTION, random_state=SEED)
    train, test = next(split.split(X, y, groups))
    predicted, std = cross_validation(X.iloc[train], y[train], groups[train])
    cv_scores = scores(y[train], predicted)
    model = fit_qspr(X.iloc[train], y[train])
    external_pred, external_std = predict_qspr(model, X.iloc[test])
    external_scores = scores(y[test], external_pred)
    # Modèle destiné aux nouvelles prédictions : réappris sur toutes les données.
    final_model = fit_qspr(X, y)
    w = final_model['weights'][final_model['selected']]
    importance = pd.DataFrame({'frequency': (w.abs() > COEF_TOL).mean(),
                               'mean_signed_weight': w.mean(),
                               'mean_abs_weight': w.abs().mean()})
    importance['relative_weight'] = importance.mean_abs_weight / importance.mean_abs_weight.sum()
    importance = importance.sort_values('relative_weight', ascending=False)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data.iloc[train].assign(predicted_cv=predicted, std_models=std).to_csv(OUTPUT/'predictions_cv.csv', index=False)
    data.iloc[test].assign(predicted_external=external_pred, std_models=external_std).to_csv(OUTPUT/'predictions_external.csv', index=False)
    pd.DataFrame([cv_scores, external_scores], index=['CV', 'external']).to_csv(OUTPUT/'metrics.csv')
    importance.to_csv(OUTPUT/'importance_l1.csv', index_label='descriptor')
    final_model['weights'].to_csv(OUTPUT/'weights_per_resample.csv', index_label='resample')
    X.to_csv(OUTPUT/'descriptors.csv', index=False)
    audit.to_csv(OUTPUT/'protonation.csv', index=False)
    plot_results(y[train], predicted, std, importance, cv_scores, OUTPUT)
    print('Validation croisée :', cv_scores)
    print('Test externe :', external_scores)
    print('UNION finale :', len(final_model['selected']), 'descripteurs')
    print('Résultats :', OUTPUT.resolve())
    return final_model


if __name__ == '__main__':
    main()
