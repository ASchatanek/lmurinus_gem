from src import PathOrganizer
from src import Model_Exploration_Tool as es

import pandas as pd

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)

lmur_pathorg = PathOrganizer()
model_dir = lmur_pathorg.get_model_version_path(mo_name="LigilactobacillusMurinus")
model = lmur_pathorg.load_model(model_dir)

kwoji = {
    "EX_cpd00001_e0": (-1000, 1000),
    "EX_cpd00009_e0": (-60, 1000),
    "EX_cpd00011_e0": (-2, 1000),
    "EX_cpd00013_e0": (-40, 1000),
    "EX_cpd00023_e0": (-2, 1000),
    "EX_cpd00028_e0": (-0.05, 1000),
    "EX_cpd00029_e0": (-50, 1000),
    "EX_cpd00030_e0": (-0.25, 1000),
    "EX_cpd00033_e0": (-1, 1000),
    "EX_cpd00035_e0": (-1, 1000),
    "EX_cpd00039_e0": (-1, 1000),
    "EX_cpd00041_e0": (-2, 1000),
    "EX_cpd00048_e0": (-76.35, 1000),
    "EX_cpd00051_e0": (-1, 1000),
    "EX_cpd00053_e0": (-2, 1000),
    "EX_cpd00054_e0": (-1, 1000),
    "EX_cpd00060_e0": (-1, 1000),
    "EX_cpd00065_e0": (-1, 1000),
    "EX_cpd00066_e0": (-1, 1000),
    "EX_cpd00067_e0": (-1000, 1000),
    "EX_cpd00069_e0": (-1, 1000),
    "EX_cpd00084_e0": (-2, 1000),
    "EX_cpd00092_e0": (-0.1, 1000),
    "EX_cpd00099_e0": (-20, 1000),
    "EX_cpd00104_e0": (-0.005, 1000),
    "EX_cpd00107_e0": (-1, 1000),
    "EX_cpd00119_e0": (-2.5, 1000),
    "EX_cpd00129_e0": (-1, 1000),
    "EX_cpd00132_e0": (-2, 1000),
    "EX_cpd00156_e0": (-1, 1000),
    "EX_cpd00161_e0": (-1, 1000),
    "EX_cpd00182_e0": (-0.1, 1000),
    "EX_cpd00184_e0": (-0.1, 1000),
    "EX_cpd00205_e0": (-90, 1000),
    "EX_cpd00218_e0": (-0.001, 1000),
    "EX_cpd00220_e0": (-0.001, 1000),
    "EX_cpd00254_e0": (-25, 1000),
    "EX_cpd00305_e0": (-0.001, 1000),
    "EX_cpd00311_e0": (-0.1, 1000),
    "EX_cpd00322_e0": (-2, 1000),
    "EX_cpd00393_e0": (-0.001, 1000),
    "EX_cpd00588_e0": (-10, 1000),
    "EX_cpd00644_e0": (-0.001, 1000),
    "EX_cpd00971_e0": (-100, 1000),
    "EX_cpd01217_e0": (-0.1, 1000),
    "EX_cpd10515_e0": (-0.2, 1000),
    "EX_cpd11606_e0": (-0.01, 1000),
    "EX_cpd00027_e0": (-1000, 1000),  # glucose
}

closed_uptake = [
    "EX_cpd00007_e0",  # Oxygen
]

essential = [
    "EX_cpd00034_e0",  # Zn2+
    "EX_cpd00058_e0",  # Cu2+
    "EX_cpd00063_e0",  # Ca2+
    "EX_cpd00149_e0",  # Co2+
    "EX_cpd10516_e0",  # fe3
]

c_sources = [
    "EX_cpd00027_e0",  # glucose
    "EX_cpd01200_e0",  # palatinose
    "EX_cpd00076_e0",  # sucrose
    "EX_cpd00082_e0",  # fructose
]

kwoji_updated = {
    "EX_cpd00588_e0": (-10, 1000),
    "EX_cpd00009_e0": (-60, 1000),
    "EX_cpd00205_e0": (-95, 1000),
    "EX_cpd00067_e0": (-92, 1000),
    "EX_cpd00254_e0": (-25, 1000),
    "EX_cpd00048_e0": (-25.56, 1000),
    "EX_cpd00001_e0": (-180.21, 1000),
    "EX_cpd00030_e0": (-0.25, 1000),
    "EX_cpd00013_e0": (-10, 1000),
    "EX_cpd00099_e0": (-27.02, 1000),
    "EX_cpd10515_e0": (-0.2, 1000),
    "EX_cpd00971_e0": (-52, 1000),
    "EX_cpd00011_e0": (-2, 1000),
    "EX_cpd00149_e0": (-0.01, 1000),
    "EX_cpd00063_e0": (-1, 1000),
    "EX_cpd00034_e0": (-0.1, 1000),
    "EX_cpd00058_e0": (-0.01, 1000),
    "EX_cpd00029_e0": (-50, 1000),
    "EX_cpd00028_e0": (-0.05, 1000),
    "EX_cpd00322_e0": (-2, 1000),
    "EX_cpd00060_e0": (-1, 1000),
    "EX_cpd00066_e0": (-1, 1000),
    "EX_cpd00119_e0": (-2.5, 1000),
    "EX_cpd00053_e0": (-2, 1000),
    "EX_cpd00023_e0": (-2, 1000),
    "EX_cpd00132_e0": (-2, 1000),
    "EX_cpd00041_e0": (-2, 1000),
    "EX_cpd00069_e0": (-1, 1000),
    "EX_cpd00065_e0": (-1, 1000),
    "EX_cpd00035_e0": (-1, 1000),
    "EX_cpd00156_e0": (-1, 1000),
    "EX_cpd00051_e0": (-1, 1000),
    "EX_cpd00129_e0": (-1, 1000),
    "EX_cpd00033_e0": (-1, 1000),
    "EX_cpd00039_e0": (-1, 1000),
    "EX_cpd00161_e0": (-1, 1000),
    "EX_cpd00107_e0": (-1, 1000),
    "EX_cpd00054_e0": (-1, 1000),
    "EX_cpd00084_e0": (-2, 1000),
    "EX_cpd00311_e0": (-0.1, 1000),
    "EX_cpd00182_e0": (-0.1, 1000),
    "EX_cpd01217_e0": (-0.1, 1000),
    "EX_cpd00184_e0": (-0.1, 1000),
    "EX_cpd00092_e0": (-0.1, 1000),
    "EX_cpd00218_e0": (-0.01, 1000),
    "EX_cpd00305_e0": (-0.01, 1000),
    "EX_cpd00644_e0": (-0.01, 1000),
    "EX_cpd00220_e0": (-0.01, 1000),
    "EX_cpd00393_e0": (-0.01, 1000),
    "EX_cpd00104_e0": (-0.005, 1000),
    "EX_cpd11606_e0": (-0.005, 1000),
    "EX_cpd00027_e0": (-1000, 1000),  # glucose
}

#     print(rxn.build_reaction_string(use_metabolite_names=False))
#     print(rxn.build_reaction_string(use_metabolite_names=True))

# for met in model.metabolites:
#     for rxn in met.reactions:
#         if "bio1" in rxn.id:
#             if len(met.reactions) == 3:
#                 print(f"{met.id} --> {met.name}")
#                 print("-" * 20)
#                 for rxn in met.reactions:
#                     if "bio1" not in rxn.id:
#                         print(f"{rxn.id}: {rxn.name}")
#                         print("Reaction:")
#                         print(rxn.build_reaction_string(use_metabolite_names=False))
#                         print(rxn.build_reaction_string(use_metabolite_names=True))
#                         print("...")
#                 print("")

exploration = es(model=model)
exploration.set_media(medium=kwoji_updated, essential=essential, closed=closed_uptake)

lala = exploration.contrained_medium_fba_analysis(percentage=0.5)

lala

test = exploration.fetch_constrained_medium_fba_fluxes(percentage=0.5)

test

from pathlib import Path

Path.cwd()
