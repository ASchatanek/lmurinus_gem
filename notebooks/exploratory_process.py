import pandas as pd

from src import PathOrganizer
from src import ExploratoryAid as EA

from notebooks import exploratory_variables as ev

################

po = PathOrganizer()
ea = EA()

LOfModels = ea.iter_through_draft_models(target_models=ev.draftModels)

LOfModels

model_stats = dict()

for model in LOfModels:

    model_stats[model.model.id] = (
        len(model.model.reactions),
        len(model.model.metabolites),
        len(model.model.genes),
    )

    print(model.model.id)
    print(f"Reactions: {len(model.model.reactions)}")
    print(f"Metabolites: {len(model.model.metabolites)}")
    print(f"Genes: {len(model.model.genes)}")
    print("")

model_stats

results = pd.DataFrame.from_dict(model_stats)

results.index = ["Reactions", "Metabolites", "Genes"]

filename = "modelsStatistics.xlsx"

results.T.sort_values(by="Reactions", ascending=False).to_excel(filename)


po = PathOrganizer()

lmurinus = po.get_model_version_path(mo_name="LigilactobacillusMurinus")

lmurinus = po.load_model(lmurinus)

lmurinus.reactions.get_by_id("rxn05629_c0")
