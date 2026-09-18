import os
from IPython.display import display, HTML, clear_output
from IPython.display import Image
import numpy as np
import pandas as pd
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


def plot_langmuir_model(model):
    time = model.root.output.solution.solution_times
    c = model.root.output.solution.unit_001.solution_outlet
    plt.plot(time, c)
    plt.title("Column (Outlet)")
    plt.xlabel("$time~/~min$")
    plt.ylabel("$concentration~/~mol \cdot L^{-1} $")
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
