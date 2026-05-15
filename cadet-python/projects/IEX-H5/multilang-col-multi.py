#!/bin/bash

import os, subprocess
from cadet import Cadet
import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv

load_dotenv()

# Cadet.cadet_path = os.getenv("cadet_root_dir")


# Define system settings

## General
init_c = [1.0, 0]
ncomp = len(init_c)  # number of components to consider in the chromatographic medium
t_cycle = 30  # in s, total time for the global unit to operate
Q = 6e-6  # in m**-3/s, flow rate
V_mixer = 0.5  # m**-3, Volume of the reactor
## Reaction
kfwd = [0.1]  # in SI
kbwd = [0.0]  # in SI
stoich = [-1, 1]  # -1 -> +1

# Initialization
model = Cadet()

# Inlet profile : https://cadet.github.io/v6.0.0a1/interface/unit_operations/inlet.html#inlet-config
## Define the inlet
model.root.input.model.unit000.unit_type = "INLET"
## Number of components in the unit (sec)
model.root.input.model.unit000.ncomp = ncomp
## Mathematic model for the intlet
model.root.input.mode.unit000.intlet_type = "PIECEISE_CUBIC_POLY"  # Only existing one
## Concentration variation in each unit in mol/m**-3
### = a + b (t-t_i) + c(t-t_i**2) + d (t-t_i**3)
model.root.input.model.unit_000.sec_000.const_coeff = [
    1.0,
]  # a
model.root.input.model.unit_000.sec_000.lin_coeff = [
    0.0,
]  # b
model.root.input.model.unit_000.sec_000.quad_coeff = [
    0.0,
]  # c
model.root.input.model.unit_000.sec_000.cube_coeff = [
    0.0,
]  # d

## Create the unit
model.root.input.model.unit_001.unit_type = "CSTR"
## number of componenets in the unit
model.root.input.model.unit_001.ncomp = ncomp
## volume f the unit
model.root.input.model.unit_001.init_volume = V_mixer
# concentration of the compenents in the unit
model.root.input.model.unit_001.init_c = [
    0.0,
    0.0,
]
## Concentration variation in each unit in mol/m**-3
### = a + b (t-t_i) + c(t-t_i**2) + d (t-t_i**3)
model.root.input.model.unit_001.sec_000.const_coeff = [
    1.0,
]  # a
model.root.input.model.unit_001.sec_000.lin_coeff = [
    0.0,
]  # b
model.root.input.model.unit_001.sec_000.quad_coeff = [
    0.0,
]  # c
model.root.input.model.unit_001.sec_000.cube_coeff = [
    0.0,
]  # d

# Columns : Build like a reactor without the reaction


# Outlet profile
## Define the outlet
model.root.input.model.unit_002.unit_type = "OUTLET"
## Number of components in the unit
model.root.input.model.unit_002.ncomp = 1
## Concentration variation in each unit in mol/m**-3
### = a + b (t-t_i) + c(t-t_i**2) + d (t-t_i**3)
model.root.input.model.unit_002.sec_001.const_coeff = [
    0.0,
]  # a
model.root.input.model.unit_002.sec_001.lin_coeff = [
    0.0,
]  # b
model.root.input.model.unit_002.sec_001.quad_coeff = [
    0.0,
]  # c
model.root.input.model.unit_002.sec_001.cube_coeff = [
    0.0,
]  # d


# Defining time passage in each section
model.root.input.solver.sections.nsec = 2  # number of sections : create discontinuit in the process (injection of specific volume/ concentration)
model.root.input.solver.sections.section_times = [
    0.0,
    t_cycle / 2,
    t_cycle,
]  # intervals between sections
model.root.input.solver.sections.section_continuity = [
    0
]  # Continuity (restart at zero between sections)
model.root.input.model.solver.user_solution_times = np.linspace(0, t_cycle, 1001)

# Connecting the element
## number of valve (connector)
model.root.input.model.connections.nswitches = 1
##
model.root.input.model.connections.switch_000.section = 0
## Matrix of connectivity is defined as : [UnitOpID from, UnitOpID to, Component from, Component to, volumetric flow rate]
model.root.input.model.connections.switch_000.connections = [
    (0, 1, -1, -1, Q),  # unit_000, unit_001, all components, all components, Q/ m^3/s
    (1, 2, -1, -1, Q),
]

# run a simulation
## save results
model.filename = "filename.h5"
## save parameters
model.save()
## run the simulation
data = model.run()
## load results
model.load()

# get the results
##
time = model.root.output.solution.solution_times
##
c = model.root.output.solution.unit_001.solution_outlet_comp_001
