import sys

sys.path.append("../")

from utils import *
# %run ../utils.py

t_cycle = 60  # 1 min = 60 s
Q = 1e-6
n_comp = 1

model = get_cadet_template(n_units=5)

model.root.input.model.unit_000.unit_type = "INLET"
model.root.input.model.unit_000.ncomp = n_comp
model.root.input.model.unit_000.inlet_type = "PIECEWISE_CUBIC_POLY"

model.root.input.model.unit_001.unit_type = "INLET"
model.root.input.model.unit_001.ncomp = n_comp
model.root.input.model.unit_001.inlet_type = "PIECEWISE_CUBIC_POLY"

model.root.input.model.unit_002.unit_type = "CSTR"
model.root.input.model.unit_002.ncomp = n_comp
model.root.input.model.unit_002.init_c = n_comp * [0.0]
model.root.input.model.unit_002.init_volume = 1e-6

model.root.input.model.unit_003.unit_type = "LUMPED_RATE_MODEL_WITH_PORES"
model.root.input.model.unit_003.ncomp = n_comp
model.root.input.model.unit_003.col_dispersion = 2.00e-7
model.root.input.model.unit_003.cross_section_area = 1.0386890710931253e-4
model.root.input.model.unit_003.col_length = 0.05
model.root.input.model.unit_003.col_porosity = 0.4
model.root.input.model.unit_003.film_diffusion = n_comp * [1e-4]


model.root.input.model.unit_00.unit_type = "OUTLET"
model.root.input.model.unit_001.ncomp = 1


model.root.input.solver.sections.nsec = 2
model.root.input.solver.sections.section_times = [0.0, t_cycle / 2, t_cycle]  # s
model.root.input.solver.sections.section_continuity = [0]
model.root.input.model.unit_000.sec_000.const_coeff = [
    1.0,
]  # mM
model.root.input.model.unit_000.sec_000.lin_coeff = [
    0.0,
]
model.root.input.model.unit_000.sec_000.quad_coeff = [
    0.0,
]
model.root.input.model.unit_000.sec_000.cube_coeff = [
    0.0,
]

model.root.input.model.unit_000.sec_001.const_coeff = [
    0.0,
]  # mM
model.root.input.model.unit_000.sec_001.lin_coeff = [
    0.0,
]
model.root.input.model.unit_000.sec_001.quad_coeff = [
    0.0,
]
model.root.input.model.unit_000.sec_001.cube_coeff = [
    0.0,
]
model.root.input.model.connections.nswitches = 1
model.root.input.model.connections.switch_000.section = 0

model.root.input.model.connections.switch_000.connections = [
    0,
    1,
    -1,
    -1,
    Q,  # unit_000, unit_001, all components, all components, Q/ m^3/s
]

model.root.input.solver.user_solution_times = np.linspace(0, t_cycle, 1001)

run_simulation(model, "test.h5")

plt.figure()

time = model.root.output.solution.solution_times
c = model.root.output.solution.unit_001.solution_inlet
plt.plot(time / 60, c)
plt.xlabel("$time~/~min$")
plt.ylabel("$Outlet~concentration~/~mM$")
plt.show()
