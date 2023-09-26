from src.data import PathOrganizer

mo_name = "Lactobacillus_murinus"

lmur_pathorg = PathOrganizer()

model_dir = lmur_pathorg.get_model_version_path()

model = lmur_pathorg.load_model(model_dir)

print(model)
print("-" * 20)
print("Model Stats")
print("-" * 20)
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
print(f"Genes: {len(model.genes)}")
print("-" * 20)

model.remove_reactions(reactions=["EX_cpd17042_e0", "EX_cpd17043_e0", "EX_cpd17041_e0"])

target = model.metabolites.get_by_id("cpd17041_e0")
target.id = "cpd17041_c0"

target = model.metabolites.get_by_id("cpd17042_e0")
target.id = "cpd17042_c0"

target = model.metabolites.get_by_id("cpd17043_e0")
target.id = "cpd17043_c0"

model.add_boundary(metabolite=model.metabolites.get_by_id("cpd17041_c0"), type="sink")
model.add_boundary(metabolite=model.metabolites.get_by_id("cpd17042_c0"), type="sink")
model.add_boundary(metabolite=model.metabolites.get_by_id("cpd17043_c0"), type="sink")

for bound in model.sinks:
    print(bound.name)
    print(bound.id)
    for met in bound.metabolites:
        print(met.id)

    print("")

# lmur_pathorg.save_model()
