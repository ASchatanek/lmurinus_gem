from src.data import PathOrganizer

mo_name = "Lactobacillus_murinus"

lmur_pathorg = PathOrganizer()

draft_path = lmur_pathorg.get_draft_model_path("lmurinus_draft.sbml")

draftmodel = lmur_pathorg.load_model(draft_path)

print(draftmodel)
print("-" * 20)
print("Model Stats")
print("-" * 20)
print(f"Reactions: {len(draftmodel.reactions)}")
print(f"Metabolites: {len(draftmodel.metabolites)}")
print(f"Genes: {len(draftmodel.genes)}")
print("-" * 20)
