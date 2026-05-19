# %% [markdown]

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

Col : Multi Langmuir and General Rate model to start and see if the first position seems realistic

Outlet with the 3 mol from Inlet2
"""

# %%
import sys

sys.path.append("..")
from utils import *

# %%
langmuir_model = get_cadet_template(n_units=5)
n_comp = 3

# col param
col_vol = 4.7 * 10**-6
col_l = 0.1
col_sec = col_vol / col_l

# DOI : 10.1002/jssc.202200943
# DOI : 10.1016/j.chroma.2004.10.008
# DOI : 10.1002/jctb.7239
# DOI : 10.1016/B978-0-323-95879-0.50003-5
col_por = 0.40
col_disp =  # Dax = u(ax) * L * sigma² / tau² and u(ax) = L / t0 (deriative of position : to time Pos/Time) or u(ax) = us / eps
par_por = 0.75
par_rad = 45e-6

#
RT = 60
t_cycle = 20 * col_vol * RT

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
langmuir_model.root.input.model.unit_002.init_volume = 1e-6
langmuir_model.root.input.model.unit_002.init_c = n_comp * [0.0]

# Column
langmuir_model.root.input.model.unit_003.unit_type = "LUMPED_RATE_MODEL_WITH_PORES"
langmuir_model.root.input.model.unit_003.ncomp = n_comp

langmuir_model.root.input.model.unit_003.col_length = col_l
langmuir_model.root.input.model.unit_003.cross_section_area = col_sec
langmuir_model.root.input.model.unit_003.col_porosity = 0.37
langmuir_model.root.input.model.unit_003.par_porosity = 0.33
langmuir_model.root.input.model.unit_003.par_radius = (90e-6)/2

langmuir_model.root.input.model.unit_003.col_dispersion = 2.0e-7
langmuir_model.root.input.model.unit_003.film_diffusion = [1e-4, 1e-4, 1e-4]

langmuir_model.root.input.model.unit_003.adsorption_model = "MULTI_COMPONENT_LANGMUIR"
langmuir_model.root.input.model.unit_003.nbound = [1] * n_comp
langmuir_model.root.input.model.unit_003.adsorption.is_kinetic = False
langmuir_model.root.input.model.unit_003.adsorption.mcl_ka = [2, 3, 1.5]
langmuir_model.root.input.model.unit_003.adsorption.mcl_kd = [2.5, 1.5, 1]
langmuir_model.root.input.model.unit_003.adsorption.mcl_qmax = [2 * 1, 1.5 * 1, 2 * 1]

langmuir_model.root.input.model.unit_003.init_c = n_comp * [0.0]
langmuir_model.root.input.model.unit_003.init_q = n_comp * [0.0]

## Outlet
langmuir_model.root.input.model.unit_004.unit_type = "OUTLET"
langmuir_model.root.input.model.unit_004.ncomp = n_comp

set_discretization(langmuir_model, n_col=40)

## Sections and Connections
langmuir_model.root.input.solver.sections.nsec = 2
langmuir_model.root.input.solver.sections.section_times = [0.0, 50.0, t_cycle]
langmuir_model.root.input.solver.sections.section_continuity = [0]

langmuir_model.root.input.model.unit_000.sec_000.const_coeff = n_comp * [1.0]
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
plt.xlabel("$time~/~min$")
plt.ylabel("$concentration~/~mM$")

time = langmuir_model.root.output.solution.solution_times
c = langmuir_model.root.output.solution.unit_002.solution_outlet
plt.figure()
plt.plot(time, c)
plt.title("Mixer (Outlet)")
plt.xlabel("$time~/~min$")
plt.ylabel("$concentration~/~mM$")

time = langmuir_model.root.output.solution.solution_times
c = langmuir_model.root.output.solution.unit_003.solution_outlet
plt.figure()
plt.plot(time, c)
plt.title("Column (Outlet)")
plt.xlabel("$time~/~min$")
plt.ylabel("$concentration~/~mM$")
