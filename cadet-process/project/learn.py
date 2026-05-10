# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: cadet_workshop
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Project template
#
# File made to run quick simulation based on a simple template easily transformable.
#
# First approach will just contain Classic chrommatography template being  :
#
# _________________________________
# ## Overall schematic for the code
#
#
# #### Mixtures → Inner Flow Circuits → Columns → Seperated molecules
#
#
# | Mixtures                             | Inner Flow Circuits  | Column     | Seperated Molecules|
# | -------------------------------------|----------------------|------------|--------------------|
# | Components                           | Flow rates           | Technology | peak 1             |
# | Fractions (molar or massic if known) | Diameters            | Ligands    | ...                |
# |                                      | Porosity of conduits | Porosity   | peak n             |
# |                                      |                      | Isotherm   |                    |
#
# __________________
#
# ## Software in use
#
# Cadet will be primarily used for the simulation, it is an open sourced software with many
# isotherms in use to run the different adsorbtion inside the chrommatographic process.
#
# It will be paired with other essential libraries such as numpy, pandas, scipy etc... Often used in science.

# %%
import numpy as np
import modin.pandas as pd # activate multi-threading for pandas
import matplotlib.pyplot as plt
import sys, os

# %%
# CADET imports

## Create process
from CADETProcess.processModel import (Inlet,Outlet,FlowSheet,Process,
                                       ComponentSystem,LumpedRateModelWithPores,
                                       MassActionLaw,
                                       Cstr, TubularReactor)
from CADETProcess.modelBuilder import SerialZone, ParallelZone


## Select a Model
"""
Linear
Langmuir
LangmuirLDF
LangmuirLDFLiquidPhase
BiLangmuir
BiLangmuirLDF
AntiLangmuir
Spreading
StericMassAction
MultistateStericMassAction
SimplifiedMultistateStericMassAction
BiStericMassAction
MobilePhaseModulator
ExtendedMobilePhaseModulator
SelfAssociation
GeneralizedIonExchange
MultiComponentColloidal
FreundlichLDF
Saska
HICWaterOnhydrophobicSurfaces
HICConstantWaterActivity
HICUnified
Sips # not yet available
"""
from CADETProcess.processModel import Langmuir

## Perform simulation
from CADETProcess.simulator import Cadet

# %%
# Indicate components lists
components_system = ComponentSystem()
components_system.add_component(["Components"])

# %% [markdown]
# [Isotherms](https://cadet.github.io/master/modelling/binding)

# %%
# Configure the model
# Size of the list must be the same as the number of components
# characteristics highly depends on the model here Langmiur is used as an example
# refer to wiki for the isotherms specifications
binding_model = Langmuir(components_system, name='langmuir')
binding_model.is_kinetic = False
binding_model.adsorption_rate = [0.02, 0.03]
binding_model.desorption_rate = [1, 1]
binding_model.capacity = [100, 100]

# %%
# Define constants if needed
Q = 1e-6/60

# %%
reaction_system = MassActionLaw(components_system)
reaction_system.add_reaction(
        indices=[0,1],
        coefficients=[-1, 1],
        k_fwd=0.1,
        k_bwd=0
        )

# %%
# Create the flowsheet

## Inner flows : Feed + Eluent
eluent = Inlet(components_system, name='eluent')
eluent.c = [0, 0]

feed = Inlet(components_system, name="feed")
feed.c = [10, 10]

## Outer flow
oulet = Outlet(components_system, name='oulet')

# %%
# Mixer
valve = Cst(components_system, 'valve')
valve.init_liquid_volume = 1e-6
valve.const_solid_volume = 0
valve.flow_rate = Q

# %%
# TubularReactor
tubing = TubularReactor(components_system, 'multi-zone')
tubing.length = 0.5
tubing.cross_section_area = 1e-5
tubing.axial_dispersion = 1e-5

# %%
# Column setup
column = LumpedRateModelWithPores(component_system, 'column')
column.binding_model = binding_model
column.length = 0.1
column.diameter = 0.01
column.bed_porosity = 0.37
column.particle_radius = 4.5e-5
column.particle_porosity = 0.33
column.axial_dispersion = 2.0e-7
column.film_diffusion = [1e-4, 1e-4]

# %%
# flow sheet
flow_sheet = FlowSheet(component_system)

flow_sheet.add_unit(feed)
flow_sheet.add_unit(eluent)
flow_sheet.add_unit(column)
flow_sheet.add_unit(outlet)
#flow_sheet.add_unit(valve)
#flow_sheet.add_unit(tubing)

flow_sheet.add_connection(feed, column)
flow_sheet.add_connection(eluent, column)
flow_sheet.add_connection(column, outlet)

"""
If a mixer is needed, (most likely we are in the real world)

# valve with tubular
flow_sheet.add_connection(feed, valve)
flow_sheet.add_connection(eluent, valve)
flow_sheet.add_connection(valve, tubing)
flow_sheet.add_connection(tubing, column)
flow_sheet.add_connection(column, outlet)
"""

# %%
# Process creation
process = Process(flow_sheet, 'batch_elution')

# events
## injecting feed,
process.add_event(name='feed_on', parameter_path='flow_sheet.feed.flow_rate', state=Q, time=0)
process.add_event('feed_off', 'flow_sheet.feed.flow_rate', 0.0, 60)

## inejecting eluent
process.add_event(name='eluent_on', parameter_path='flow_sheet.eluent.flow_rate', state=Q)
process.add_event('eluent_off', 'flow_sheet.eluent.flow_rate', 0.0)

# dependancies
process.add_event_dependency('eluent_on', 'feed_off')
process.add_event_dependency('eluent_off', 'feed_on')

## Set process times
process.cycle_time = 1200 # in seconds

# %%
simulation_results = process_simulator.simulate(process)
_ = simulation_results.solution.column.outlet.plot()



