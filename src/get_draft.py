from src.data import PathOrganizer

mo_name = "Lactobacillus_murinus"

lmur_pathorg = PathOrganizer()

test = lmur_pathorg.get_draft_model_path("lmurinus_draft.sbml")

draftmodel = lmur_pathorg.load_model(test)

print(draftmodel)
