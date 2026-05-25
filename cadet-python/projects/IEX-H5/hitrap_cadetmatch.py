# %%
import sys
from dotenv import load_dotenv

sys.path.append("..")
from utils import *
from CADETProcess.reference import ReferenceIO
from CADETProcess.comparison import Comparator
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import FlowSheet
from CADETProcess.processModel import Process
from CADETProcess.processModel import Langmuir
from CADETProcess.processModel import Inlet, Cstr, LumpedRateModelWithPores, Outlet
from CADETProcess.simulator import Cadet
from CADETProcess.optimization import OptimizationProblem

load_dotenv()

# %%

# akta config
akta_mixer_v = 1.4 * 10**-6

# col param
col_vol = 5 * 10**-6  # m**3
col_l = 0.025  # m
col_sec = col_vol / col_l  # m**2
par_por = 0.75  # find the ref
par_rad = 45e-6
RT = 60  # s
t_cycle = 37 * col_vol * 10**6 * RT

# experimental results
Q = (5 * 10**-6) / 60  # m**3/s
dead_time = float(os.getenv("IEXH5_HITRAP_T0")) * 60  # s
dead_vol = float(os.getenv("IEXH5_HITRAP_V0")) * 10**-6  # m**3
sigma = float(os.getenv("IEXH5_HITRAP_sigma"))
tau = float(os.getenv("IEXH5_HITRAP_tau")) * 60  # s
data = f"{os.getenv('data')}/IEX-HYDRO/Hitrap_Capto_Q_for_HA_29-01-2026_pH_8_2_001.xlsx"

W = 4 * sigma
N = 16 * (tau / W) ** 2
H = col_l / N
u = col_l / dead_time

# DOI : 10.1002/jssc.202200943
# DOI : 10.1016/j.chroma.2004.10.008
# DOI : 10.1002/jctb.7239
# DOI : 10.1016/B978-0-323-95879-0.50003-5
# DOI : carta's book
# DOI : 10.1063/5.0276752
# DOI : 10.1016/S0021-9673(03)00311-X

v = col_l / dead_time
u = Q / col_sec
col_por = u / v
print(col_por)

# %%
uax = col_l / dead_time
col_disp = H * u / 2

# protein information
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
mw_impu2 = 3.8  # kda
rm_impu2 = compute_rm(mw_impu2)

mw = [mw_prot1, mw_impu1, mw_impu2]
rm = [rm_prot1, rm_impu1, rm_impu2]

# theorical

data = f"{os.getenv('data')}/IEX-HYDRO/Hitrap_Capto_Q_for_HA_29-01-2026_pH_8_2_001.xlsx"
hitrap_data = pd.read_excel(data)
x = [i / 5 for i in hitrap_data["ml"]]
y = hitrap_data["mAU"]

reference = ReferenceIO("ref for captoQ", x, y)
reference.plot()
comparator = Comparator()
comparator.add_reference(reference)
comparator.add_difference_metric("RMSE", reference, "column.outlet")

Q = (5**10 - 6) / 60

# Component System
component_system = ComponentSystem(["protein", "impurity1", "impurity2"])

# Binding Model
binding_model = Langmuir(component_system, name="langmuir")
binding_model.is_kinetic = False
binding_model.adsorption_rate = [0.02, 0.03, 0.05]
binding_model.desorption_rate = [1, 1, 1]
binding_model.capacity = [100, 100, 100]

# Unit Operations
feed = Inlet(component_system, name="feed")
feed.c = [1, 1, 1]

eluent = Inlet(component_system, name="eluent")
eluent.c = [0, 0, 0]

# Mixer
valve = Cstr(component_system, "valve")
valve.init_liquid_volume = 1.4**10 - 6
valve.const_solid_volume = 0
valve.flow_rate = Q

# Column
column = LumpedRateModelWithPores(component_system, "column")
column.binding_model = binding_model
column.length = col_l
column.diameter = col_sec
column.bed_porosity = col_por
column.particle_radius = par_rad
column.particle_porosity = par_por
column.axial_dispersion = col_disp
column.film_diffusion = [1e-4, 1e-4, 1e-4]

outlet = Outlet(component_system, name="outlet")

# Flow Sheet
flow_sheet = FlowSheet(component_system)

flow_sheet.add_unit(feed)
flow_sheet.add_unit(eluent)
flow_sheet.add_unit(valve)
flow_sheet.add_unit(column)
flow_sheet.add_unit(outlet)

flow_sheet.add_connection(feed, valve)
flow_sheet.add_connection(eluent, valve)
flow_sheet.add_connection(valve, column)
flow_sheet.add_connection(column, outlet)

# Process
process = Process(flow_sheet, "batch elution")

## Create Events and Durations
process.add_event("feed_on", "flow_sheet.feed.flow_rate", Q, 0)
process.add_event("feed_off", "flow_sheet.feed.flow_rate", 0.0, t_cycle / 2)

process.add_event("eluent_off", "flow_sheet.eluent.flow_rate", 0.0, 0.0)
process.add_event("eluent_on", "flow_sheet.eluent.flow_rate", Q, 60)

## Set Process Times
process.cycle_time = t_cycle

process_simulator = Cadet()
simulation_results = process_simulator.simulate(process)
_ = simulation_results.solution.column.inlet.plot()
_ = simulation_results.solution.column.outlet.plot()

# comparator.plot_comparison(simulation_results)

optimization_problem = OptimizationProblem("porosity_axial_dispersion")

optimization_problem.add_evaluation_object(process)
optimization_problem.add_variable(
    name="ka",
    parameter_path="flow_sheet.binding_model.adsorption_rate",
    lb=0.1,
    ub=0.8,
    transform="auto",
)
optimization_problem.add_variable(
    name="axial_dispersion",
    parameter_path="flow_sheet.column.axial_dispersion",
    lb=1e-9,
    ub=1e-6,
    transform="auto",
)
optimization_problem.add_variable(
    name="ka",
    parameter_path="flow_sheet.binding_model.desorption_rate",
    lb=0.1,
    ub=0.8,
    transform="auto",
)
optimization_problem.add_variable(
    name="ka",
    parameter_path="flow_sheet.binding_model.capacity",
    lb=0.1,
    ub=0.8,
    transform="auto",
)
optimization_problem.add_variable(
    name="ka",
    parameter_path="flow_sheet.column.film_diffusion",
    lb=0.1,
    ub=0.8,
    transform="auto",
)

simulator = Cadet()

optimization_problem.add_evaluator(simulator)

optimization_problem.add_objective(
    comparator, n_objectives=comparator.n_metrics, requires=[simulator]
)
