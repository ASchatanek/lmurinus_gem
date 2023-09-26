from src.data import PathOrganizer
import re

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

for bound in model.sinks:
    print(bound.name)
    print(bound.id)
    for met in bound.metabolites:
        print(met.id)

    print("")

# lmur_pathorg.save_model()

print(len(model.boundary))
print(len(model.exchanges))
print(len(model.sinks))
print(type(model.demands))

for exchange in model.exchanges:
    ex_reac = exchange.reaction
    met_id = re.match(r"([a-z]+)(\d+)\w([a-z]\d)", ex_reac).group()
    met = model.metabolites.get_by_id(met_id)

    match = re.match(r"(\w+)-e0-e0", met.name)
    if match:
        print(f"{met} --> {met.name}")

        met.name = f"{match.group(1)}-e0"

        print(f"{met} --> {met.name}")


model.metabolites.get_by_id("cpd00158_e0")

lmur_pathorg.save_model()
