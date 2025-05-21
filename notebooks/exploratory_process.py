import set_cwd
import notebooks.exploratory_variables as ev


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

draftModels = [
    "amuciniphila_draft_xml.xml",
    "amucosicola_draft_xml.xml",
    "bcaecimuris_draft_xml.xml",
    "bpseudococcoides_draft_xml.xml",
    "cinnocuum_draft_xml.xml",
    "cramosum_draft_xml.xml",
    "eclostridioformis_draft_xml.xml",
    "efaecalis_draft_xml.xml",
    "ecoli_draft_xml.xml",
    "emuris_draft_xml.xml",
    "fplautii_draft_xml.xml",
    "fbutyricus_draft_xml.xml",
    "lreuteri_draft_xml.xml",
    "lmurinus_draft_xml.xml",
    "mschaedleri_draft_xml.xml",
    "mintestinale_draft_xml.xml",
    "pgoldsteinii_draft_xml.xml",
    "tmuris_draft_xml.xml",
    "xrodentium_draft_xml.xml",
]

##############################

for model in draftModels:

    m_dir = po.get_draft_model_path(
        draft_name=model,
    )

    draftModel = po.load_model(m_dir)

    es_dModel = es(
        model=draftModel,
    )

    print(draftModel)

    print("---")

    essential, variable = es_dModel.find_crucial_and_variable_essentialMetabolites(
        baseMedium=ev.kwoji_medium,
        baseMediumDataframe=ev.kwoji_medium_df,
        closedMetabolites=ev.closed_uptake,
    )

    print("Essential:")

    for ess in essential:

        print(draftModel.reactions.get_by_id(ess).name)

    es_dModel.set_media(
        medium=ev.kwoji_medium,
        closed=ev.closed_uptake,
        essential=essential,
    )

    print(draftModel.slim_optimize())

    print("---")

    for rxn in variable:

        yeah = essential + [rxn]

        es_dModel.set_media(
            medium=ev.kwoji_medium,
            closed=ev.closed_uptake,
            essential=yeah,
        )

        print(
            f"{draftModel.reactions.get_by_id(rxn).name} --> {draftModel.slim_optimize()}"
        )

    print("_______________________________")

##############################

test = EA()

gotit = test.metabolite_assay(
    target_models=ev.problematicModels,
    medium=ev.kwoji_medium,
    closed=ev.closed_uptake,
    essential=ev.essential,
    metabolites_df=ev.biolog_met_df,
    medium_df=ev.kwoji_medium_df,
)

gotit

###############################
am_dir = po.get_draft_model_path(
    draft_name="ecoli_draft_xml.xml",
)
am_m = po.load_model(am_dir)

ex_am = es(model=am_m)

ex_am.add_and_set_media(
    medium=ev.kwoji_medium,
    closed_m=ev.closed_uptake,
    essential_m=None,
    medium_df=ev.kwoji_medium_df,
)

am_m.slim_optimize()
am_m.reactions.get_by_id("EX_cpd00025_e0")

essential, possibly = ex_am.find_essentialMetabolites(
    baseMedium=ev.kwoji_medium,
    baseMediumDataframe=ev.kwoji_medium_df,
    closedMetabolites=ev.closed_uptake,
)

essential

possibly

am_m.reactions.get_by_id("EX_cpd00025_e0")

ex_am.add_and_set_media(
    medium=ev.kwoji_medium,
    closed_m=ev.closed_uptake,
    essential_m=essential,
    medium_df=ev.kwoji_medium_df,
)

print(am_m.slim_optimize())

for p in possibly:

    esemero = essential + [p]

    ex_am.add_and_set_media(
        medium=ev.kwoji_medium,
        closed_m=ev.closed_uptake,
        essential_m=esemero,
        medium_df=ev.kwoji_medium_df,
    )
    print(am_m.reactions.get_by_id(p).name)
    print(am_m.slim_optimize())
    print("-------------")

######################

ex_am.add_and_set_media(
    medium=ev.kwoji_medium,
    closed_m=ev.closed_uptake,
    essential_m=ev.essential,
    medium_df=ev.kwoji_medium_df,
)

am_m.summary().uptake_flux.sort_values(["flux"], ascending=False)

max_growth = am_m.slim_optimize()
testingMinMedia = cobra.medium.minimal_medium(
    model=am_m,
    min_objective_value=max_growth,
    open_exchanges=False,
)

print(max_growth)
for x in testingMinMedia.index:
    rxnName = am_m.reactions.get_by_id(x)
    min_flux = testingMinMedia[x]

    print(f"{rxnName.name} --> {min_flux}")

#####################


def create_ATP_maintenance_reaction():

    reaction = cobra.Reaction(
        id="rxn11300_c0",
        name="ATP maintenance",
        lower_bound=0.0,
        upper_bound=1000.0,
    )

    reaction.add_metabolites(
        {
            am_m.metabolites.get_by_id("cpd00002_c0"): -1.0,
            am_m.metabolites.get_by_id("cpd00001_c0"): -1.0,
            am_m.metabolites.get_by_id("cpd00008_c0"): 1.0,
            am_m.metabolites.get_by_id("cpd00009_c0"): 1.0,
            am_m.metabolites.get_by_id("cpd00067_c0"): 1.0,
        }
    )

    return reaction


with am_m:

    reaction = create_ATP_maintenance_reaction()

    am_m.add_reactions([reaction])

    am_m.objective = reaction

    print(cobra.flux_analysis.pfba(am_m).objective_value)
    print(am_m.optimize().objective_value)

am_m.optimize().objective_value

am_m.metabolites.cpd00002_c0.summary()

am_m.metabolites.get_by_id("cpd00067_c0")

am_m.reactions.get_by_id("rxn00216_c0")

for rxn in am_m.reactions:

    print(f"{rxn.id} --- {rxn.name}")

###############################
## LMurinus
lm_dir = po.get_draft_model_path(draft_name="lmurinus_draft_xml.xml")
lm_m = po.load_model(lm_dir)
ex_lm = es(model=lm_m)

lm_m

ex_lm.add_missing_medium_met(medium=ev.kwoji_medium)

for m in lm_m.metabolites:

    print(m)

ex_lm.set_media(
    medium=ev.kwoji_medium,
    essential=ev.essential,
    closed=ev.closed_uptake,
)

r_lm = ex_lm.add_assay_metabolites(mets_df=ev.biolog_met_df)

r_lm

###############################
## Pgoldsteinii
pg_dir = po.get_draft_model_path(draft_name="pgoldsteinii_draft_xml.xml")
pgold_m = po.load_model(pg_dir)
ex_pgold = es(model=pgold_m)

ex_pgold.add_missing_medium_met(medium=ev.kwoji_medium)

ex_pgold.set_media(
    medium=ev.kwoji_medium,
    essential=ev.essential,
    closed=ev.closed_uptake,
)

results = ex_pgold.add_assay_metabolites(mets_df=ev.biolog_met_df)

results

###############################

test = EA()

test.metabolite_assay(
    target_models=ev.draftModels,
    medium=ev.kwoji_medium,
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

ex_pgold.find_medium_outliers(up, medium=ev.kwoji_medium)
ex_pgold.find_medium_outliers(sec, medium=ev.kwoji_medium)

# Lmurinus
model_dir = po.get_model_version_path(mo_name="LigilactobacillusMurinus")
model = po.load_model(model_dir)

exploration = es(model=model)
exploration.set_media(
    medium=ev.kwoji_medium,
    essential=ev.essential,
    closed=ev.closed_uptake,
)

model.summary()

exploration.display_hm_constrained_medium_fba_analysis(
    percentage=0.5,
    cmap="Spectral",
)
