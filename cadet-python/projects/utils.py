import os
from IPython.display import display, HTML, clear_output
from IPython.display import Image
import numpy as np
import modin.pandas as pd
import scipy
from addict import Dict
import json
import matplotlib.pyplot as plt
from ipywidgets import interact, interactive
import ipywidgets as widgets
import tempfile
import subprocess
from CADETMatch.jupyter import Match
import shutil
import platform
from pathlib import Path
from cadet import Cadet


def get_cadet_template(n_units=3, split_components_data=False):
    cadet_template = Cadet()

    cadet_template.root.input.model.nunits = n_units

    # Store solution
    cadet_template.root.input["return"].split_components_data = split_components_data
    cadet_template.root.input["return"].split_ports_data = 0
    cadet_template.root.input["return"].unit_000.write_solution_inlet = 1
    cadet_template.root.input["return"].unit_000.write_solution_outlet = 1
    cadet_template.root.input["return"].unit_000.write_solution_bulk = 1
    cadet_template.root.input["return"].unit_000.write_solution_particle = 1
    cadet_template.root.input["return"].unit_000.write_solution_solid = 1
    cadet_template.root.input["return"].unit_000.write_solution_flux = 1
    cadet_template.root.input["return"].unit_000.write_solution_volume = 1
    cadet_template.root.input["return"].unit_000.write_coordinates = 1
    cadet_template.root.input["return"].unit_000.write_sens_outlet = 1

    for unit in range(n_units):
        cadet_template.root.input["return"]["unit_{0:03d}".format(unit)] = (
            cadet_template.root.input["return"].unit_000
        )

    # Tolerances for the time integrator
    cadet_template.root.input.solver.time_integrator.abstol = 1e-6
    cadet_template.root.input.solver.time_integrator.algtol = 1e-10
    cadet_template.root.input.solver.time_integrator.reltol = 1e-6
    cadet_template.root.input.solver.time_integrator.init_step_size = 1e-6
    cadet_template.root.input.solver.time_integrator.max_steps = 1000000

    # Solver settings
    cadet_template.root.input.model.solver.gs_type = 1
    cadet_template.root.input.model.solver.max_krylov = 0
    cadet_template.root.input.model.solver.max_restarts = 10
    cadet_template.root.input.model.solver.schur_safety = 1e-8

    # Run the simulation on single thread
    cadet_template.root.input.solver.nthreads = 1

    return cadet_template


def set_discretization(model, n_col=20, n_par_types=1):
    columns = {
        "GENERAL_RATE_MODEL",
        "LUMPED_RATE_MODEL_WITH_PORES",
        "LUMPED_RATE_MODEL_WITHOUT_PORES",
    }

    for unit_name, unit in model.root.input.model.items():
        if "unit_" in unit_name and unit.unit_type in columns:
            unit.discretization.spatial_method = "FV"
            unit.discretization.ncol = n_col
            unit.discretization.npar = 5
            unit.discretization.npartype = n_par_types

            unit.discretization.par_disc_type = "EQUIDISTANT_PAR"
            unit.discretization.use_analytic_jacobian = 1
            unit.discretization.reconstruction = "WENO"
            unit.discretization.gs_type = 1
            unit.discretization.max_krylov = 0
            unit.discretization.max_restarts = 10
            unit.discretization.schur_safety = 1.0e-8

            unit.discretization.weno.boundary_model = 0
            unit.discretization.weno.weno_eps = 1e-10
            unit.discretization.weno.weno_order = 3


def create_langmuir_model():

    langmuir_model = get_cadet_template(n_units=3)

    # INLET
    langmuir_model.root.input.model.unit_000.unit_type = "INLET"
    langmuir_model.root.input.model.unit_000.ncomp = 1
    langmuir_model.root.input.model.unit_000.inlet_type = "PIECEWISE_CUBIC_POLY"

    # Column
    langmuir_model.root.input.model.unit_001.unit_type = "LUMPED_RATE_MODEL_WITH_PORES"
    langmuir_model.root.input.model.unit_001.ncomp = 1

    langmuir_model.root.input.model.unit_001.col_length = 0.25
    langmuir_model.root.input.model.unit_001.cross_section_area = 1.0386890710931253e-4
    langmuir_model.root.input.model.unit_001.col_porosity = 0.37
    langmuir_model.root.input.model.unit_001.par_porosity = 0.33
    langmuir_model.root.input.model.unit_001.par_radius = 4.5e-5

    langmuir_model.root.input.model.unit_001.col_dispersion = 2.0e-7
    langmuir_model.root.input.model.unit_001.film_diffusion = [
        2.0e-7,
    ]

    langmuir_model.root.input.model.unit_001.adsorption_model = (
        "MULTI_COMPONENT_LANGMUIR"
    )

    langmuir_model.root.input.model.unit_001.nbound = [
        1,
    ]
    langmuir_model.root.input.model.unit_001.adsorption.is_kinetic = 1
    langmuir_model.root.input.model.unit_001.adsorption.mcl_ka = [
        1.144,
    ]
    langmuir_model.root.input.model.unit_001.adsorption.mcl_kd = [
        2.0e-3,
    ]
    langmuir_model.root.input.model.unit_001.adsorption.mcl_qmax = [
        4.88,
    ]

    langmuir_model.root.input.model.unit_001.init_c = [
        0.0,
    ]
    langmuir_model.root.input.model.unit_001.init_q = [
        0.0,
    ]

    ## Outlet
    langmuir_model.root.input.model.unit_002.ncomp = 1
    langmuir_model.root.input.model.unit_002.unit_type = "OUTLET"

    set_discretization(langmuir_model)

    # Sections and connections
    langmuir_model.root.input.solver.sections.nsec = 2
    langmuir_model.root.input.solver.sections.section_times = [0.0, 3000.0, 9500.0]
    langmuir_model.root.input.solver.sections.section_continuity = [
        0,
    ]

    ## Inlet Profile
    langmuir_model.root.input.model.unit_000.sec_000.const_coeff = [
        1.0,
    ]
    langmuir_model.root.input.model.unit_000.sec_001.const_coeff = [
        0.0,
    ]

    ## Switches
    langmuir_model.root.input.model.connections.nswitches = 1
    langmuir_model.root.input.model.connections.switch_000.section = 0
    langmuir_model.root.input.model.connections.switch_000.connections = [
        0,
        1,
        -1,
        -1,
        2.88e-8,
        1,
        2,
        -1,
        -1,
        2.88e-8,
    ]

    # set the times that the simulator writes out data for
    langmuir_model.root.input.solver.user_solution_times = np.linspace(0, 9500, 9501)

    return langmuir_model


def plot_langmuir_model(model):
    time = model.root.output.solution.solution_times
    c = model.root.output.solution.unit_001.solution_outlet
    plt.plot(time, c)
    plt.title("Column (Outlet)")
    plt.xlabel("$time~/~min$")
    plt.ylabel("$concentration~/~mol \cdot L^{-1} $")
    plt.show()


def create_dextran_model():

    dextran_model = get_cadet_template(n_units=3, split_components_data=True)

    # INLET
    dextran_model.root.input.model.unit_000.unit_type = "INLET"
    dextran_model.root.input.model.unit_000.ncomp = 1
    dextran_model.root.input.model.unit_000.inlet_type = "PIECEWISE_CUBIC_POLY"

    # Column
    dextran_model.root.input.model.unit_001.unit_type = "LUMPED_RATE_MODEL_WITH_PORES"
    dextran_model.root.input.model.unit_001.ncomp = 1
    dextran_model.root.input.model.unit_001.nbound = [0]

    dextran_model.root.input.model.unit_001.col_length = 0.25
    dextran_model.root.input.model.unit_001.cross_section_area = 1e-4
    dextran_model.root.input.model.unit_001.col_porosity = 0.37
    dextran_model.root.input.model.unit_001.par_porosity = 0.33
    dextran_model.root.input.model.unit_001.par_radius = 4.5e-5

    dextran_model.root.input.model.unit_001.col_dispersion = 2.0e-7
    dextran_model.root.input.model.unit_001.film_diffusion = [
        0.0,
    ]

    dextran_model.root.input.model.unit_001.adsorption_model = "NONE"

    dextran_model.root.input.model.unit_001.init_c = [
        0.0,
    ]

    set_discretization(dextran_model, n_col=100)

    ## Outlet
    dextran_model.root.input.model.unit_002.ncomp = 1
    dextran_model.root.input.model.unit_002.unit_type = "OUTLET"

    # Sections and connections
    dextran_model.root.input.solver.sections.nsec = 2
    dextran_model.root.input.solver.sections.section_times = [0.0, 50.0, 600.0]
    dextran_model.root.input.solver.sections.section_continuity = [
        0,
    ]

    ## Inlet Profile
    dextran_model.root.input.model.unit_000.sec_000.const_coeff = [
        1.0,
    ]
    dextran_model.root.input.model.unit_000.sec_001.const_coeff = [
        0.0,
    ]

    ## Switches
    dextran_model.root.input.model.connections.nswitches = 1
    dextran_model.root.input.model.connections.switch_000.section = 0
    dextran_model.root.input.model.connections.switch_000.connections = [
        0,
        1,
        -1,
        -1,
        2.88e-8,
        1,
        2,
        -1,
        -1,
        2.88e-8,
    ]

    # set the times that the simulator writes out data for
    dextran_model.root.input.solver.user_solution_times = np.linspace(0, 600, 601)

    return dextran_model


def plot_dextran_model(model):
    time = model.root.output.solution.solution_times
    c = model.root.output.solution.unit_001.solution_outlet_comp_000
    plt.plot(time, c)
    plt.title("Column (Outlet)")
    plt.xlabel("$time~/~min$")
    plt.ylabel("$concentration~/~mol$")
    plt.show()


def run_simulation(cadet, file_name=None):
    cadet.filename = file_name
    # save the simulation
    cadet.save()

    # run the simulation and load results
    data = cadet.run()
    cadet.load()

    # Raise error if simulation fails
    if data.return_code == 0:
        print("Simulation completed successfully")
    else:
        print(data)
        raise Exception("Simulation failed")
