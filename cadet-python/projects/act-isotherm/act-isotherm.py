# %%
# imports
import os
from cadet import Cadet
import numpy as np
import matplotlib.pyplot as plt

# %%
home_directory = os.environ.get("HOME")
cadet_folder = "Work/cadet"
cadet_bin = "CADET-Core/install/bin/cadet-cli"
cadet_normal = "cadet-python"
cadet_bin_path = home_directory + "/" + cadet_folder + "/" + cadet_bin
savefig_path = (
    home_directory + "/" + cadet_folder + "/" + cadet_normal + "/act-isotherm/graphs"
)

print(
    f"HOME DIR {home_directory} \n",
    f"CADET FOLDER {cadet_folder} \n",
    f"CADET BIN {cadet_bin} \n",
    f"CADET NORMAL {cadet_normal} \n",
    f"CADET BIN PATH {cadet_bin_path} \n",
    f"SAVEFIG PATH {savefig_path} \n",
)
Cadet.cadet_path = f"{cadet_bin_path}"

# %%
col_order = [
    "5_mL_MSS",
    "5_mL_CAP",
    "5_mL_HC650F",
    "1_mL_CAP",
    "1_mL_MSS",
    "5_mL_MSS",
    "5_mL_CAP",
    "5_mL_HC650F",
]

# %%
tables = [
    ## Table1
    {
        "name": col_order[0],
        "mAb": "B",
        0.392: {
            "residence_time_min": [0.84, 0.91, 1.43, 2.50, 1.43, 1.43, 1.43],
            "gradient_length_CV": [15, 12, 10, 7, 5, 20, 30],
            "pH_wash_pH_elute": [
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
            ],
        },
        0.731: {
            "residence_time_min": [0.84, 0.91, 1.12, 2.50, None, None, None],
            "gradient_length_CV": [15, 12, 10, 7, None, None, None],
            "pH_wash_pH_elute": [
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                None,
                None,
                None,
            ],
        },
        1.216: {
            "residence_time_min": [0.84, 0.91, 1.12, 2.50, None, None, None],
            "gradient_length_CV": [15, 12, 10, 7, None, None, None],
            "pH_wash_pH_elute": [
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                None,
                None,
                None,
            ],
        },
    },
    ## Table2
    {
        "name": col_order[1],
        "mAb": "B",
        0.240: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
            ],
        },
        0.746: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
            ],
        },
        2.334: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
            ],
        },
    },
    ## Table3
    {
        "name": col_order[2],
        "mAb": "B",
        0.312: {
            "residence_time_min": [2.00, 2.50, 3.00, 3.85, 5.00, 6.00, 6.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [5.2, 3.0],
                [5.2, 3.0],
                [5.2, 3.0],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
            ],
        },
        0.558: {
            "residence_time_min": [2.00, 2.50, 3.00, 3.85, 5.00, 6.00, 6.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [5.2, 3.0],
                [5.2, 3.0],
                [5.2, 3.0],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
            ],
        },
        1.498: {
            "residence_time_min": [2.00, 2.50, 3.00, 3.85, 5.00, 6.00, 6.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [5.2, 3.0],
                [5.2, 3.0],
                [5.2, 3.0],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
            ],
        },
    },
    ## Table4
    {
        "name": col_order[3],
        "mAb": "B",
        0.097: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
            ],
        },
        0.298: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
            ],
        },
        0.749: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
            ],
        },
    },
    ## Table5
    {
        "name": col_order[4],
        "mAb": "B",
        0.097: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
            ],
        },
        0.298: {
            "residence_time_min": [1.00, 1.43, 2.00, 3.00, 4.00, 5.00, 5.00],
            "gradient_length_CV": [20, 15, 12, 10, 10, 7, 5],
            "pH_wash_pH_elute": [
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
            ],
        },
    },
    ## Table6
    {
        "name": col_order[5],
        "mAb": "A",
        0.250: {
            "residence_time_min": [1.43, 1.00, 2.00, 3.00, 4.00, 5.00, 3.33, 2.00],
            "gradient_length_CV": [30, 20, 15, 10, 5, 0, 0, 10],
            "pH_wash_pH_elute": [
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
            ],
        },
        0.620: {
            "residence_time_min": [1.43, 1.00, 2.00, 3.00, 4.00, 5.00, 3.33, 2.00],
            "gradient_length_CV": [30, 20, 15, 10, 5, 0, 0, 10],
            "pH_wash_pH_elute": [
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
            ],
        },
        2.521: {
            "residence_time_min": [1.43, 1.00, 2.00, 3.00, 4.00, 5.00, None, None],
            "gradient_length_CV": [30, 20, 15, 10, 5, 0, None, None],
            "pH_wash_pH_elute": [
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.0],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
                [None, None],
                [None, None],
            ],
        },
    },
    ## Table7
    {
        "name": col_order[6],
        "mAb": "A",
        0.265: {
            "residence_time_min": [1.43, 1.00, 2.00, 3.00, 4.00, 5.00],
            "gradient_length_CV": [30, 20, 15, 10, 5, 0],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
            ],
        },
        2.695: {
            "residence_time_min": [1.43, 1.00, 2.00, 3.00, 4.00, 5.00],
            "gradient_length_CV": [30, 20, 15, 10, 5, 0],
            "pH_wash_pH_elute": [
                [6.0, 3.0],
                [6.0, 3.0],
                [6.0, 3.0],
                [5.5, 3.3],
                [5.5, 3.3],
                [5.5, 3.3],
            ],
        },
    },
    # Table8
    {
        "name": col_order[7],
        "mAb": "A",
        0.265: {
            "residence_time_min": [2.50, 2.00, 3.33, 3.85, 5.00, 6.00],
            "gradient_length_CV": [30, 20, 15, 10, 5, 0],
            "pH_wash_pH_elute": [
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
            ],
        },
        0.743: {
            "residence_time_min": [2.50, 2.00, 3.33, 3.85, 5.00, 6.00],
            "gradient_length_CV": [30, 20, 15, 10, 5, 0],
            "pH_wash_pH_elute": [
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
                [5.2, 2.7],
            ],
        },
    },
]


# %%
## System settings
def create_column(name: str, mAb: str) -> tuple:
    global col_order

    """
    Quickly define column parameters
    """

    if name == col_order[0]:
        column_length = 0.025  ## m
        column_radius = 0.016  ## m
        column_volume = ((np.pi * column_radius**2) / 4) * column_length  ## m^3
        column_porosity = 0.31  ## 1
        particle_porosity = 0.96  ## 1
        total_porosity = column_porosity + (1.0 - column_porosity) * particle_porosity
        if mAb == "A":
            qmax = 55.6  ## mg/mL
            keq = 102  ## mL/mg
            dp_e5 = 1.99
            aeta = 1.81
            etaG = 2.28
            pKaA = 2.07
            pKaG = 5.29
        elif mAb == "B":
            qmax = 55.6
            keq = 102
            dp_e5 = 1.99
            aeta = 1.81
            etaG = 2.28
            pKaA = 3.82
            pKaG = 5.29
        else:
            raise ValueError(f"mAb asked not available {mAb}")

    elif name == col_order[1]:
        column_length = 0.050
        column_radius = 0.0113
        column_volume = ((np.pi * column_radius**2) / 4) * column_length
        column_porosity = 0.37
        particle_porosity = 0.97
        total_porosity = column_porosity + (1.0 - column_porosity) * particle_porosity
        if mAb == "A":
            qmax = 55.6
            keq = 102
            dp_e5 = 1.99
            aeta = 1.81
            etaG = 2.28
            pKaA = 2.07
            pKaG = 5.29
        elif mAb == "B":
            qmax = 50.3
            keq = 268
            dp_e5 = 1.25
            aeta = 2.12
            etaG = 3.48
            pKaA = 3.82
            pKaG = 4.66
        else:
            raise ValueError("mAb asked not available")

    elif name == col_order[2]:
        column_length = 0.100
        column_radius = 0.008
        column_volume = ((np.pi * column_radius**2) / 4) * column_length
        column_porosity = 0.34
        particle_porosity = 0.91
        total_porosity = column_porosity + (1.0 - column_porosity) * particle_porosity
        if mAb == "A":
            qmax = 80.2
            keq = 320
            dp_e5 = 1.62
            aeta = 3.50
            etaG = 2.41
            pKaA = 3.13
            pKaG = 5.11
        elif mAb == "B":
            qmax = 80.7
            keq = 542
            dp_e5 = 1.65
            aeta = 1.55
            etaG = 4.64
            pKaA = 1.05
            pKaG = 4.33
        else:
            raise ValueError("mAb asked not available")

    elif name == col_order[3]:
        column_length = 0.020
        column_radius = 0.008
        column_volume = ((np.pi * column_radius**2) / 4) * column_length
        column_porosity = 0.37
        particle_porosity = 0.97
        total_porosity = column_porosity + (1.0 - column_porosity) * particle_porosity
        if mAb == "A":
            qmax = 46.4
            keq = 115
            dp_e5 = None
            aeta = None
            etaG = None
            pKaA = None
            pKaG = None
        elif mAb == "B":
            qmax = 51
            keq = 134
            dp_e5 = 0.88
            aeta = 2.48
            etaG = 2.48
            pKaA = 4.10
            pKaG = 4.10
        else:
            raise ValueError("mAb asked not available")

    elif name == col_order[4]:
        column_length = 0.025
        column_radius = 0.007
        column_volume = ((np.pi * column_radius**2) / 4) * column_length
        column_porosity = 0.31
        particle_porosity = 0.96
        total_porosity = column_porosity + (1.0 - column_porosity) * particle_porosity
        if mAb == "A":
            qmax = 48.1
            keq = 301
            dp_e5 = None
            aeta = None
            etaG = None
            pKaA = None
            pKaG = None
        elif mAb == "B":
            qmax = 50.5
            keq = 291
            dp_e5 = 0.72
            aeta = 0.82
            etaG = 3.33
            pKaA = 5.93
            pKaG = 4.34
        else:
            raise ValueError("mAb asked not available")

    else:
        raise ValueError(f"Name is not defined : {name}")
    return (
        column_length,
        column_volume,
        column_porosity,
        particle_porosity,
        total_porosity,
        qmax,
        keq,
        dp_e5,
        aeta,
        etaG,
        pKaA,
        pKaG,
    )


# %%
for i, table in enumerate(tables):
    mab = table.get("mAb")
    for loading in list(
        table.keys()
    ):  # loading will be used as is in the parameter definition of the col (packing concentration)
        if loading == "mAb":
            continue

        if loading == "name":
            (
                column_length,
                column_volume,
                column_porosity,
                particle_porosity,
                total_porosity,
                qmax,
                keq,
                dp_e5,
                n_proton_a,
                n_proton_g,
                aciditea,
                aciditeg,
            ) = create_column(table[loading], mab)
            continue

        residence_time_min = table[loading]["residence_time_min"]
        gradient_length_CV = table[loading]["gradient_length_CV"]
        pH_wash_pH_elute = table[loading]["pH_wash_pH_elute"]
        groups = len(residence_time_min)

        for group in range(groups):
            if (
                residence_time_min[group] is None
                or gradient_length_CV[group] is None
                or pH_wash_pH_elute[group] is None
                or None in pH_wash_pH_elute[group]
            ):
                print(f"Skipping invalid group {group}")
                continue

            protein_MW = 150  ## kDa
            c_feed = 4.0 / protein_MW  ## mol/m^3

            RT = residence_time_min[group] * 60
            Q = column_volume / RT

            elution_pH_start = pH_wash_pH_elute[group][0]
            elution_pH_end = pH_wash_pH_elute[group][1]

            G = gradient_length_CV[group]
            load_cv = G * 0.2
            wash_cv = G * 0.3
            elution_cv = G * 0.5
            event_CV = [
                0,
                load_cv,
                wash_cv,
                elution_cv,
            ]  ## load, wash, elution in CV

            print(
                f"protein MW: {protein_MW} \n",
                f"Q: {Q} \n",
                f"c_feed: {c_feed} \n",
                f"RT: {RT} \n",
                f"elution_pH_start: {elution_pH_start} \n",
                f"elution_pH_end: {elution_pH_end} \n",
                f"event_CV: {event_CV} \n",
                f"column_length: {column_length} \n",
                f"column_volume: {column_volume} \n",
                f"column_porosity: {column_porosity} \n",
                f"particle_porosity: {particle_porosity} \n",
                f"total_porosity: {total_porosity} \n",
                f"qmax: {qmax} \n",
                f"keq: {keq} \n",
                f"dp_e5: {dp_e5} \n",
                f"n_proton_a: {n_proton_a} \n",
                f"n_proton_g: {n_proton_g} \n",
                f"aciditea: {aciditea} \n",
                f"aciditeg: {aciditeg} \n",
            )

            # %%
            model = Cadet()

            model.root.input.model.nunits = 3

            model.root.input.model.unit_000.unit_type = "INLET"
            model.root.input.model.unit_000.ncomp = 2  ## The first component is pH
            model.root.input.model.unit_000.inlet_type = "PIECEWISE_CUBIC_POLY"

            model.root.input.model.unit_001.unit_type = "GENERAL_RATE_MODEL"
            model.root.input.model.unit_001.ncomp = 2

            ## Geometry
            model.root.input.model.unit_001.col_length = column_length  # m
            model.root.input.model.unit_001.cross_section_area = (
                column_volume / column_length
            )  # m^2
            model.root.input.model.unit_001.col_porosity = column_porosity  # 1
            model.root.input.model.unit_001.par_porosity = particle_porosity  # 1
            model.root.input.model.unit_001.par_radius = 0.0425e-3  # m

            ## Transport
            model.root.input.model.unit_001.col_dispersion = 1.36e-8  # m^2/s
            model.root.input.model.unit_001.film_diffusion = [
                1,
                1.41e-5,
            ]  # m/s
            model.root.input.model.unit_001.par_diffusion = [1, 1.99e-11]  # m^2/s
            model.root.input.model.unit_001.par_surfdiffusion = [0.0, 0.0]

            model.root.input.model.unit_001.adsorption_model = (
                "AFFINITY_COMPLEX_TITRATION"
            )

            model.root.input.model.unit_001.adsorption.is_kinetic = 1
            model.root.input.model.unit_001.adsorption.act_ka = [
                1.0,
                keq * protein_MW * (1.0 - total_porosity),
            ]  ##  m^3 solid phase / mol protein / s   ml/mg=m^3/kg    1 kda = 1 kg/mol
            model.root.input.model.unit_001.adsorption.act_kd = [1.0, 1.0]  ## s^-1
            model.root.input.model.unit_001.adsorption.act_qmax = [
                1e-10,
                qmax / protein_MW / (1.0 - total_porosity),
            ]  ##  mol/m^3 solid phase mg/ml=kg/m^3 / 150 kg/mol = 1/150 mol/m^3

            model.root.input.model.unit_001.adsorption.act_etaA = [
                0,
                n_proton_a,
            ]
            model.root.input.model.unit_001.adsorption.act_pkaA = [
                0,
                aciditea,
            ]
            model.root.input.model.unit_001.adsorption.act_etaG = [
                0,
                n_proton_g,
            ]
            model.root.input.model.unit_001.adsorption.act_pkaG = [
                0,
                aciditeg,
            ]

            model.root.input.model.unit_001.init_c = [
                7.2,
                0.0,
            ]
            model.root.input.model.unit_001.init_q = [0.0, 0.0]

            ### Grid cells
            model.root.input.model.unit_001.discretization.ncol = 50
            model.root.input.model.unit_001.discretization.npar = 12

            ### Bound states
            model.root.input.model.unit_001.discretization.nbound = [
                0,
                1,
            ]

            ### Other options
            model.root.input.model.unit_001.discretization.par_disc_type = (
                "EQUIDISTANT_PAR"
            )
            model.root.input.model.unit_001.discretization.use_analytic_jacobian = 1
            model.root.input.model.unit_001.discretization.reconstruction = "WENO"
            model.root.input.model.unit_001.discretization.gs_type = 1
            model.root.input.model.unit_001.discretization.max_krylov = 0
            model.root.input.model.unit_001.discretization.max_restarts = 10
            model.root.input.model.unit_001.discretization.schur_safety = 1.0e-8

            model.root.input.model.unit_001.discretization.weno.boundary_model = 0
            model.root.input.model.unit_001.discretization.weno.weno_eps = 1e-10
            model.root.input.model.unit_001.discretization.weno.weno_order = 2

            model.root.input.model.unit_002.unit_type = "OUTLET"
            model.root.input.model.unit_002.ncomp = 2

            model.root.input.solver.sections.nsec = 4
            model.root.input.solver.sections.section_times = [
                0,
                event_CV[1] * RT,
                (event_CV[1] + event_CV[2]) * RT,
                (event_CV[1] + event_CV[2] + event_CV[3]) * RT,
                (event_CV[1] + event_CV[2] + event_CV[3] + 5.5) * RT,
            ]  ## s
            model.root.input.solver.sections.section_continuity = [0, 0, 0, 0]

            ## load
            model.root.input.model.unit_000.sec_000.const_coeff = [
                7.2,
                c_feed,
            ]  # mol / m^3
            model.root.input.model.unit_000.sec_000.lin_coeff = [
                0.0,
                0.0,
            ]
            model.root.input.model.unit_000.sec_000.quad_coeff = [
                0.0,
                0.0,
            ]
            model.root.input.model.unit_000.sec_000.cube_coeff = [
                0.0,
                0.0,
            ]

            ## wash
            model.root.input.model.unit_000.sec_001.const_coeff = [
                5.5,
                0.0,
            ]
            model.root.input.model.unit_000.sec_001.lin_coeff = [
                0.0,
                0.0,
            ]
            model.root.input.model.unit_000.sec_001.quad_coeff = [
                0.0,
                0.0,
            ]
            model.root.input.model.unit_000.sec_001.cube_coeff = [
                0.0,
                0.0,
            ]

            ## elution
            model.root.input.model.unit_000.sec_002.const_coeff = [
                elution_pH_start,
                0.0,
            ]
            model.root.input.model.unit_000.sec_002.lin_coeff = [
                -(elution_pH_start - elution_pH_end) / (event_CV[3] * RT),
                0.0,
            ]
            model.root.input.model.unit_000.sec_002.quad_coeff = [
                0.0,
                0.0,
            ]
            model.root.input.model.unit_000.sec_002.cube_coeff = [
                0.0,
                0.0,
            ]

            ## end
            model.root.input.model.unit_000.sec_003.const_coeff = [elution_pH_end, 0.0]
            model.root.input.model.unit_000.sec_003.lin_coeff = [
                0.0,
                0.0,
            ]
            model.root.input.model.unit_000.sec_003.quad_coeff = [
                0.0,
                0.0,
            ]
            model.root.input.model.unit_000.sec_003.cube_coeff = [
                0.0,
                0.0,
            ]

            model.root.input.model.connections.nswitches = 1
            model.root.input.model.connections.switch_000.section = 0
            model.root.input.model.connections.switch_000.connections = [
                0,
                1,
                -1,
                -1,
                Q,
                1,
                2,
                -1,
                -1,
                Q,
            ]

            model.root.input.model.solver.gs_type = 1
            model.root.input.model.solver.max_krylov = 0
            model.root.input.model.solver.max_restarts = 10
            model.root.input.model.solver.schur_safety = 1e-8

            model.root.input.solver.nthreads = 1

            # Tolerances for the time integrator
            model.root.input.solver.time_integrator.abstol = 1e-6
            model.root.input.solver.time_integrator.algtol = 1e-10
            model.root.input.solver.time_integrator.reltol = 1e-6
            model.root.input.solver.time_integrator.init_step_size = 1e-6
            model.root.input.solver.time_integrator.max_steps = 1000000

            # Return data
            model.root.input["return"].split_components_data = 1
            model.root.input["return"].split_ports_data = 0
            model.root.input["return"].unit_000.write_solution_bulk = 0
            model.root.input["return"].unit_000.write_solution_inlet = 0
            model.root.input["return"].unit_000.write_solution_outlet = 1

            # Copy settings to the other unit operations
            model.root.input["return"].unit_001 = model.root.input["return"].unit_000
            model.root.input["return"].unit_002 = model.root.input["return"].unit_000

            # Solution times
            model.root.input.solver.user_solution_times = np.linspace(
                0, (event_CV[1] + event_CV[2] + event_CV[3] + 5.5) * RT, 201
            )

            model.filename = f"h5/act-isotherms-{table.get('mAb')}-{col_order[i]}.h5"
            model.save()
            data = model.run()

            model.load()

            time = model.root.output.solution.solution_times
            print(time, Q)
            c = model.root.output.solution.unit_001.solution_outlet_comp_001
            pH_outlet = model.root.output.solution.unit_001.solution_outlet_comp_000

            fig = plt.figure()
            title = f"Simulation using ACT for pro A using mAb {table.get('mAb')} and {col_order[i]}"
            ax = fig.add_subplot()
            ax.plot(
                time[1:] / 60 * Q * 6e7, c[1:] * protein_MW, c="orange", label="mAb"
            )
            ax.set_xlabel(r"Volume/mL")
            ax.set_ylabel(r"Concentration/(mg/mL)")

            ax_ph = ax.twinx()
            ax_ph.plot(time / 60 * Q * 6e7, pH_outlet, label="pH")
            ax_ph.set_ylabel("pH")
            plt.savefig(f"{savefig_path}/{table.get('mAb')}_{col_order[i]}.png")
# %%
