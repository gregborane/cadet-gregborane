import os, subprocess
from cadet import Cadet
import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv

load_dotenv()

# Cadet.cadet_path = os.getenv("cadet_root_dir")

# %%
# Define system settings

## General
column_length = 0.025  ## m
column_volume = 5.025e-6  ## m^3
protein_MW = 150  ## kDa
column_porosity = 0.31  ## 1
particle_porosity = 0.95  ## 1
Q = 1.25 / (6e7)  ## m**-3/s
elution_pH_start = 5.5
elution_pH_end = 3.3
loading = 0.097  ## mg/mL
c_feed = 4.0 / protein_MW  ## mol/m^3
RT = column_volume / Q  ## s
event_CV = [
    0.0,
    loading / (c_feed * protein_MW),
    4.0,
    10.0,
]  ## load, wash, elution in CV
total_porosity = column_porosity + (1.0 - column_porosity) * particle_porosity
ncomp = 3
V_mixer = 10000000
t_cycle = 45 * 60
# %%

# Initialization
model = Cadet()

# Create Inlet
model.root.input.model.unit000.unit_type = "INLET"
model.root.input.model.unit000.ncomp = ncomp
model.root.input.mode.unit000.intlet_type = "PIECEISE_CUBIC_POLY"  # Only existing one
model.root.input.model.unit_000.sec_000.const_coeff = [
    1.0,
]
model.root.input.model.unit_000.sec_000.lin_coeff = [
    0.0,
]
model.root.input.model.unit_000.sec_000.quad_coeff = [
    0.0,
]
model.root.input.model.unit_000.sec_000.cube_coeff = [
    0.0,
]

## Create dumper
model.root.input.model.unit_001.unit_type = "CSTR"
model.root.input.model.unit_001.ncomp = ncomp
model.root.input.model.unit_001.init_volume = V_mixer
model.root.input.model.unit_001.init_c = [
    0.0,
    0.0,
]
model.root.input.model.unit_001.sec_000.const_coeff = [
    1.0,
]
model.root.input.model.unit_001.sec_000.lin_coeff = [
    0.0,
]
model.root.input.model.unit_001.sec_000.quad_coeff = [
    0.0,
]
model.root.input.model.unit_001.sec_000.cube_coeff = [
    0.0,
]

# Create the column and its properties
model.root.inut.model.unit_002.unit_type = "GENERAL_RATE_MODEL"
model.root.input.model.unit_002.ncomp = ncomp
model.root.input.model.unit_002.cross_section_area = column_volume / column_length
model.root.input.model.unit_002.col_length = column_length
model.root.input.model.unit_002.col_porosity = column_porosity  # 1
model.root.input.model.unit_002.col_dispersion = 1.36e-8  # m^2/s
model.root.input.model.unit_002.par_porosity = particle_porosity  # 1
model.root.input.model.unit_002.par_radius = 0.0425e-3  # m
model.root.input.model.unit_002.par_diffusion = [1, 1, 1]
model.root.input.model.unit_002.par_surfdiffusion = [0.0, 0.0]
model.root.input.model.unit_002.film_diffusion = [
    1,
    1.41e-5,
]
model.root.input.model.unit_002.init_c = [9, 9, 9]
model.root.input.model.unit_002.init_q = [9, 9, 9]

# Create resolution grid
model.root.input.model.unit_001.discretization.ncol = 50
model.root.input.model.unit_001.discretization.npar = 12
model.root.input.model.unit_001.discretization.nbound = [
    0,
    1,
]
model.root.input.model.unit_001.discretization.par_disc_type = "EQUIDISTANT_PAR"
model.root.input.model.unit_001.discretization.use_analytic_jacobian = 1
model.root.input.model.unit_001.discretization.reconstruction = "WENO"
model.root.input.model.unit_001.discretization.gs_type = 1
model.root.input.model.unit_001.discretization.max_krylov = 0
model.root.input.model.unit_001.discretization.max_restarts = 10
model.root.input.model.unit_001.discretization.schur_safety = 1.0e-8
model.root.input.model.unit_001.discretization.weno.boundary_model = 0
model.root.input.model.unit_001.discretization.weno.weno_eps = 1e-10
model.root.input.model.unit_001.discretization.weno.weno_order = 2

# Create Isotherm
model.root.input.model.unit_002.adsorption_model = "MULTI_COMPONENT_LANGMUIR"
model.root.input.model.unit_002.adsorption_model_multiplex = 1
model.root.input.model.unit_002.is_kinetic = 1
model.root.input.model.unit_002.mcl_ka = [0, 0, 0]
model.root.input.model.unit_002.mcl_kd = [0, 0, 0]

# Outlet profile
model.root.input.model.unit_002.unit_type = "OUTLET"
model.root.input.model.unit_002.ncomp = [1, 0.0]
model.root.input.model.unit_002.sec_001.lin_coeff = [
    0.0,
]
model.root.input.model.unit_002.sec_001.quad_coeff = [
    0.0,
]
model.root.input.model.unit_002.sec_001.cube_coeff = [
    0.0,
]

# Defining time passage in each section
model.root.input.solver.sections.nsec = 2
model.root.input.solver.sections.section_times = [
    0.0,
    t_cycle / 2,  # wash
    t_cycle,  # elution
]
model.root.input.solver.sections.section_continuity = [0]
model.root.input.model.solver.user_solution_times = np.linspace(0, t_cycle, 1001)

# Connecting the element
model.root.input.model.connections.nswitches = 1
model.root.input.model.connections.switch_000.section = 0
model.root.input.model.connections.switch_000.connections = [
    (0, 1, -1, -1, Q),
    (1, 2, -1, -1, Q),
    (2, 3, -1, -1, Q),
    (3, 4, -1, -1, Q),
]

# Run a simulation
model.filename = "simulation_of_.h5"
model.save()
data = model.run()
model.load()
time = model.root.output.solution.solution_times
c = model.root.output.solution.unit_001.solution_outlet_comp_001

# %%
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(time[1:] / 60 * Q * 6e7, c[1:] * protein_MW, c="orange", label="mAb")
ax.set_xlabel(r"Volume/mL")
ax.set_ylabel(r"Concentration/(mg/mL)")
ax_ph = ax.twinx()
ax_ph.set_ylabel("pH")
plt.show()
