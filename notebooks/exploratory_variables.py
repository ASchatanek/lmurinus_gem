from src import Reader
import pandas as pd


pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)

##############

biolog_met_df = Reader().find_dataframe(
    file_name="Metabolites Listing.xlsx",
    sheet_name="BIOLOG",
)

biolog_met_df = biolog_met_df.loc[
    :,
    ("Metabolite", "ModelSeed", "Formula"),
]

##############

kwoji_medium_df = Reader().find_dataframe(
    file_name="Metabolites Listing.xlsx",
    sheet_name="KWOJI",
)

kwoji_medium_df = kwoji_medium_df.loc[
    :,
    ("Metabolite", "ModelSeed", "Formula"),
]

##############

# * Media Definition

## * Kwoji Media
closed_uptake = [
    "EX_cpd00007_e0",  # Oxygen
]

essential = [
    "EX_cpd10516_e0",  # fe3
    "EX_cpd00264_e0",  # spermidine
    "EX_cpd00232_e0",  # Neu5Ac
    "EX_cpd00246_e0",  # Inosine
    "EX_cpd00242_e0",  # Carbonate
    "EX_cpd00076_e0",  # * Sucrose
    "EX_cpd00108_e0",  # * Galactose
    "EX_cpd00122_e0",  # N-Acetyl-D-glucosamine
    "EX_cpd00215_e0",  # Pyridoxal
    "EX_cpd00020_e0",  # Pyruvate
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

kwoji_medium = {
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
}

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

problematicModels = [
    "amuciniphila_draft_xml.xml",
    "amucosicola_draft_xml.xml",
    "bcaecimuris_draft_xml.xml",
    "bpseudococcoides_draft_xml.xml",
    "cinnocuum_draft_xml.xml",
    "cramosum_draft_xml.xml",
    "eclostridioformis_draft_xml.xml",
    "efaecalis_draft_xml.xml",
]
