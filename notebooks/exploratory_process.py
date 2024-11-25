import set_cwd
import pandas_settings
from src import PathOrganizer
from src import Model_Exploration_Tool as es
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
%matplotlib inline

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)

lmur_pathorg = PathOrganizer()
model_dir = lmur_pathorg.get_model_version_path(mo_name="LigilactobacillusMurinus")
model = lmur_pathorg.load_model(model_dir)

closed_uptake = [
    "EX_cpd00007_e0",  # Oxygen
]

essential = [
    "EX_cpd10516_e0",  # fe3
]

c_sources = [
    "EX_cpd00027_e0",  # glucose
    "EX_cpd01200_e0",  # palatinose
    "EX_cpd00138_e0",  # mannose
    "EX_cpd00076_e0",  # sucrose
    "EX_cpd00082_e0",  # fructose
    "EX_cpd00020_e0",  # pyruvate
    "EX_cpd00751_e0",  # fucose
    "EX_cpd00280_e0",  # galacturonic acid
    "EX_cpd00020_e0",  # pyruvate
]

kwoji_updated = {
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
    "EX_cpd00027_e0": (-50, 1000), # glucose
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

exploration.display_hm_constrained_medium_fba_analysis(percentage=0.5, cmap="Spectral")

model.summary()

# model.summary()

# lalista = [["EX_cpd00122_e0"],
# ["EX_cpd00366_e0"],
# ["EX_cpd00158_e0"],
# ["EX_cpd00082_e0"], 
# ["EX_cpd00751_e0"],
# ["EX_cpd00108_e0"],
# ["EX_cpd00280_e0"], 
# ["EX_cpd00222_e0"], 
# ["EX_cpd00027_e0"], 
# ["EX_cpd00100_e0"], 
# ["EX_cpd00080_e0"], 
# ["EX_cpd00121_e0"], 
# ["EX_cpd00208_e0"],
# ["EX_cpd00179_e0"],
# ["EX_cpd01262_e0"], 
# ["EX_cpd00314_e0"], 
# ["EX_cpd00138_e0"], 
# ["EX_cpd03198_e0"], 
# ["EX_cpd01200_e0"], 
# ["EX_cpd00382_e0"], 
# ["EX_cpd00396_e0"], 
# ["EX_cpd01030_e0"], 
# ["EX_cpd00588_e0"],
# ["EX_cpd00076_e0"], 
# ["EX_cpd00794_e0"],
# ["EX_cpd00029_e0"],
# ["EX_cpd00047_e0"],
# ["EX_cpd00106_e0"],
# ["EX_cpd00380_e0"],
# ["EX_cpd00221_e0"],
# ["EX_cpd00159_e0"],
# ["EX_cpd00130_e0"],
# ["EX_cpd00141_e0"],
# ["EX_cpd00020_e0"],
# ["EX_cpd00036_e0"],
# ["EX_cpd00035_e0"],
# ["EX_cpd00132_e0"],
# ["EX_cpd00023_e0"],
# ["EX_cpd00053_e0"],
# ["EX_cpd00060_e0"],
# ["EX_cpd00066_e0"],
# ["EX_cpd00054_e0"],
# ["EX_cpd00161_e0"],
# ["EX_cpd00156_e0"],
# ["EX_cpd00041_e0"],
# ["EX_cpd00246_e0"],
# ["EX_cpd00184_e0"]]

# model.metabolites.get_by_id("cpd00492_c0")
# model.reactions.get_by_id("rxn00897_c0")

# exploration.print_reactions_from_metabolite(target_metabolite_id="cpd03194_c0")

for rxn in model.reactions:
    print(f"{rxn.id} -> {rxn.name}")
    
for metabolite in model.metabolites:
    print(f"{metabolite.id} -> {metabolite.name}")


model.metabolites.get_by_id("cpd00150_e0")

model.reactions.get_by_id("rxn09979_c0")

exploration.print_reactions_from_metabolite(target_metabolite_id="cpd00150_c0")

# model.reactions.get_by_id("rxn15138_c0")
# model.optimize()
# original = model.summary().to_frame()
# original = original[original["flux"]!=0].index

# exploration.set_media(medium=kwoji_updated, essential=essential, closed=closed_uptake)
# complete_result = dict()
# for source in c_sources:
    
#     with model:
#         target = model.reactions.get_by_id(source)
#         target.bounds = (-50, 1000)
        
#         model.optimize()
#         result = model.summary().to_frame()
#         result = result[result["flux"] != 0]
        
#         print(model.optimize().objective_value)
#         res_ind = result.index
        
#         solution = result.loc[:,["flux"]].squeeze().rename(source)
        
#         complete_result[source] = solution

# final_df = pd.DataFrame(complete_result)

# nono = final_df[~final_df.index.isin(original)]
# nono = nono.dropna(how="all",axis=1)

# # model.reactions.get_by_id("EX_cpd00076_e0")

# names = dict()
# for i in nono.index: 
    
#     target = model.reactions.get_by_id(i)
    
#     for m in target.metabolites:
#         names[i] = m.name
        
# nono = nono.rename(index=names)

# co_names = dict()
# for i in nono.columns: 
    
#     target = model.reactions.get_by_id(i)
    
#     for m in target.metabolites:
#         co_names[i] = m.name

# nono = nono.rename(columns=co_names)

# nono