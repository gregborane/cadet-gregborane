# %% [markdwn]

"""

number are to be verified.
```
-sec_000------ Inlet1 ---|
                         |
-sec_001------ Inlet2 ----->---sec_001--[ Col ]-sec_002--[ Outlet ]

```

Inlet1 contains :

buffer + salt  : NaPO4(-) + [Na+;Cl-] at
pH = 8 and
[phosphate] = 50 mM
[NaCl] = 200 mM

Inlet2 contains :
Protein :
Impurities1 :
Impurities2 :

Dumper : simulate reality but might be removed

Col : Multi Langmuir and Lumped Rate model with pores to start and see if the first position seems realistic

Outlet with the 3 mol from Inlet2
"""

# %%
import sys
from dotenv import load_dotenv

sys.path.append("..")
from utils import *

load_dotenv()


# %%

# akta config
akta_mixer_v = 0.6 * 10**-6

# col param
col_vol = 5 * 10**-6  # m**3
col_l = 0.025  # m
col_sec = col_vol / col_l  # m**2
par_por = 0.75  # find the ref
par_rad = 45e-6

# experimental results
Q = (5 * 10**-6) / 60  # m**3/s
dead_time = float(os.getenv("IEXH5_HITRAP_T0")) * 60  # s
dead_vol = float(os.getenv("IEXH5_HITRAP_V0")) * 10**-6  # m**3
sigma = float(os.getenv("IEXH5_HITRAP_sigma"))
tau = float(os.getenv("IEXH5_HITRAP_tau")) * 60  # s

# %%

# DOI : 10.1002/jssc.202200943
# DOI : 10.1016/j.chroma.2004.10.008
# DOI : 10.1002/jctb.7239
# DOI : 10.1016/B978-0-323-95879-0.50003-5
# DOI : carta's book
# DOI : 10.1063/5.0276752

col_por = col_vol / dead_vol  # porosity = Vcol/DeadVol
uax = col_l / dead_time
col_disp = (
    (uax * col_l * sigma**2) / (2 * tau**2)
)  # Dax = u(ax) * L * sigma² / tau² and u(ax) = L / t0 (deriative of position : to time Pos/Time) or u(ax) = us / eps

# %%

rho_prot = 1.29 * 1000
avogadro_numb = 6.023 * 10**24


# is it concentrated ?
# is it denatured ?
# eta * Do = kb * T / (6 * pi * rH)
#
def compute_kf(rm, u, D0, Dm):
    global eta, rho, por_rad
    """
    DOI : 10.1016/j.cej.2011.07.035
    """
    ScTerm = 1.1 * ((eta / rho) * D0) ** (1 / 3)
    ReynoldTerm = ((rho * u * 2 * rm) / eta) ** 0.6
    kf = 2 + ScTerm * ReynoldTerm * (Dm / por_rad)
    return kf


def compute_rm(Mw):
    global rho_prot, avogadro_numb

    return ((3 / (4 * np.pi)) * (Mw / (rho_prot * avogadro_numb))) ** (1 / 3)


mw_prot1 = 316  # kda
rm_prot1 = compute_rm(mw_prot1)
mw_impu1 = 105  # kda
rm_impu1 = compute_rm(mw_impu1)
mw_impu2 = 3.8  # kda / projected value, not in range of calibration
rm_impu2 = compute_rm(mw_impu2)

mw = [mw_prot1, mw_impu1, mw_impu2]
rm = [rm_prot1, rm_impu1, rm_impu2]

RT = 60  # s
t_cycle = 20 * col_vol * 10**6 * RT

# %%
langmuir_model = get_cadet_template(n_units=5)
n_comp = 3

# %%

# Feed
langmuir_model.root.input.model.unit_000.unit_type = "INLET"
langmuir_model.root.input.model.unit_000.ncomp = n_comp
langmuir_model.root.input.model.unit_000.inlet_type = "PIECEWISE_CUBIC_POLY"

# INLET
langmuir_model.root.input.model.unit_001.unit_type = "INLET"
langmuir_model.root.input.model.unit_001.ncomp = n_comp
langmuir_model.root.input.model.unit_001.inlet_type = "PIECEWISE_CUBIC_POLY"

# Mixer
langmuir_model.root.input.model.unit_002.unit_type = "CSTR"
langmuir_model.root.input.model.unit_002.ncomp = n_comp
langmuir_model.root.input.model.unit_002.init_volume = akta_mixer_v
langmuir_model.root.input.model.unit_002.init_c = n_comp * [0.0]

# Column
langmuir_model.root.input.model.unit_003.unit_type = "LUMPED_RATE_MODEL_WITH_PORES"
langmuir_model.root.input.model.unit_003.ncomp = n_comp

langmuir_model.root.input.model.unit_003.col_length = col_l
langmuir_model.root.input.model.unit_003.cross_section_area = col_sec
langmuir_model.root.input.model.unit_003.col_porosity = col_por
langmuir_model.root.input.model.unit_003.par_porosity = par_por
langmuir_model.root.input.model.unit_003.par_radius = par_rad

langmuir_model.root.input.model.unit_003.col_dispersion = col_disp
langmuir_model.root.input.model.unit_003.film_diffusion = [1e-4, 1e-4, 1e-4]

langmuir_model.root.input.model.unit_003.adsorption_model = "MULTI_COMPONENT_LANGMUIR"
langmuir_model.root.input.model.unit_003.nbound = [1, 1, 1, 0, 0]
langmuir_model.root.input.model.unit_003.adsorption.is_kinetic = False
langmuir_model.root.input.model.unit_003.adsorption.mcl_ka = [
    1000 * 10**-3,
    1000 * 0.5 * 10**-3,
    1000 * 0.3 * 10**-3,
    0,
    0,
]
langmuir_model.root.input.model.unit_003.adsorption.mcl_kd = [
    1000 * 10**-3,
    1000 * 0.5 * 10**-3,
    1000 * 0.3 * 10**-3,
    0,
    0,
]
langmuir_model.root.input.model.unit_003.adsorption.mcl_qmax = [
    1000 * 10**-3,
    1000 * 4 * 10**-3,
    1000 * 0.2 * 10**-3,
    0,
    0,
]

langmuir_model.root.input.model.unit_003.init_c = n_comp * [0.0]
langmuir_model.root.input.model.unit_003.init_q = n_comp * [0.0]

## Outlet
langmuir_model.root.input.model.unit_004.unit_type = "OUTLET"
langmuir_model.root.input.model.unit_004.ncomp = n_comp

set_discretization(langmuir_model, n_col=40)

## Sections and Connections
langmuir_model.root.input.solver.sections.nsec = 2
langmuir_model.root.input.solver.sections.section_times = [0.0, 20 * 60, 37 * 60]
langmuir_model.root.input.solver.sections.section_continuity = [0]

langmuir_model.root.input.model.unit_000.sec_000.lin_coeff = [/]
langmuir_model.root.input.model.unit_001.sec_000.const_coeff = n_comp * [0.0]

langmuir_model.root.input.model.connections.nswitches = 2
langmuir_model.root.input.model.connections.switch_000.section = 0
langmuir_model.root.input.model.connections.switch_000.connections = [
    (0, 2, -1, -1, col_vol / RT),
    (2, 3, -1, -1, col_vol / RT),
    (3, 4, -1, -1, col_vol / RT),
]

langmuir_model.root.input.model.connections.switch_001.section = 1
langmuir_model.root.input.model.connections.switch_001.connections = [
    (1, 2, -1, -1, col_vol / RT),
    (2, 3, -1, -1, col_vol / RT),
    (3, 4, -1, -1, col_vol / RT),
]

# set the times that the simulator writes out data for
langmuir_model.root.input.solver.user_solution_times = np.linspace(0, t_cycle, 601)

# %%
langmuir_model.filename = f"protein_imp_simu.h5"
langmuir_model.save()
data = langmuir_model.run()
langmuir_model.load()

time = langmuir_model.root.output.solution.solution_times
c = langmuir_model.root.output.solution.unit_002.solution_inlet
plt.figure()
plt.plot(time, c)
plt.title("Mixer (Inlet)")
plt.xlabel("$time~/~sec$")
plt.ylabel("$concentration~/~mM$")

time = langmuir_model.root.output.solution.solution_times
c = langmuir_model.root.output.solution.unit_002.solution_outlet
plt.figure()
plt.plot(time, c)
plt.title("Mixer (Outlet)")
plt.xlabel("$time~/~sec$")
plt.ylabel("$concentration~/~mM$")

time = langmuir_model.root.output.solution.solution_times
c = langmuir_model.root.output.solution.unit_003.solution_outlet
plt.figure()
plt.plot(time, c)
plt.title("Column (Outlet)")
plt.xlabel("$time~/~sec$")
plt.ylabel("$concentration~/~mM$")
