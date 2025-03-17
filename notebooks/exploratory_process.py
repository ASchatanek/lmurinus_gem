import exploratory_variables as ev

# import pandas_settings
from src import PathOrganizer
from src import Model_Exploration_Tool as es
from src import ExploratoryAid as EA

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import cobra
import re

# %matplotlib inline

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)

po = PathOrganizer()

##############################

test = EA()

test.metabolite_assay(
    target_models=ev.draftModels,
    medium=ev.kwoji_updated,
    closed=ev.closed_uptake,
    essential=ev.essential,
    metabolites=ev.biolog_met_df,
)

###############################

## LMurinus
lm_dir = po.get_draft_model_path(draft_name="lmurinus_draft_xml.xml")
lm_m = po.load_model(lm_dir)
ex_lm = es(model=lm_m)

ex_lm.add_missing_medium_met(medium=ev.kwoji_updated)

ex_lm.set_media(
    medium=ev.kwoji_updated,
    essential=ev.essential,
    closed=ev.closed_uptake,
)

lm_m.summary()

lm_m.metabolites.get_by_id("cpd00246_e0")

r_lm = ex_lm.add_assay_metabolites(mets_df=ev.biolog_met_df)

r_lm

## Pgoldsteinii
pg_dir = po.get_draft_model_path(draft_name="pgoldsteinii_draft_xml.xml")
pgold_m = po.load_model(pg_dir)
ex_pgold = es(model=pgold_m)

ex_pgold.add_missing_medium_met(medium=ev.kwoji_updated)

ex_pgold.set_media(
    medium=ev.kwoji_updated,
    essential=ev.essential,
    closed=ev.closed_uptake,
)

results = ex_pgold.add_assay_metabolites(mets_df=ev.biolog_met_df)

results

###############

test = EA()

test.metabolite_assay(
    target_models=ev.draftModels,
    medium=ev.kwoji_updated,
    closed=ev.closed_uptake,
    essential=ev.essential,
    metabolites=ev.biolog_met_df,
)


for met in pgold_m.reactions.get_by_id(id="bio1").metabolites:

    print(met.id)
    print(met.name)
    print("")

pgold_m.metabolites.get_by_id(id="cpd00076_c0")
pgold_m.reactions.get_by_id(id="rxn09658_c0")

ex_pgold.print_reactions_from_metabolite(target_metabolite_id="cpd00027_e0")


up, sec = ex_pgold.gather_media_fluxes()

ex_pgold.find_medium_outliers(up, medium=ev.kwoji_updated)
ex_pgold.find_medium_outliers(sec, medium=ev.kwoji_updated)

# Lmurinus
model_dir = po.get_model_version_path(mo_name="LigilactobacillusMurinus")
model = po.load_model(model_dir)

exploration = es(model=model)
exploration.set_media(
    medium=ev.kwoji_updated,
    essential=ev.essential,
    closed=ev.closed_uptake,
)

model.summary()

exploration.display_hm_constrained_medium_fba_analysis(
    percentage=0.5,
    cmap="Spectral",
)
