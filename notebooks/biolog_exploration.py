import set_cwd
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

pd.set_option("display.precision", 3)
pd.set_option("display.max_columns", 400)
pd.set_option("display.max_rows", 400)

from src import Reader, Cleaner, BEAT, BNT, BAA

#######################

reader = Reader()
biolog_values = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="BIOLOG")
abs_535 = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="Abs_535")
abs_575 = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="Abs_575")
abs_590 = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="Abs_590")
abs_750 = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="Abs_750")

biolog_compounds_df = reader.find_dataframe(
    "plate_readout_numbers.xlsx", sheet_name="Compounds"
)

biolog_compounds = reader.get_compounds_list(biolog_compounds_df)

biolog_ids = reader.find_dataframe(
    "plate_readout_categorizations.xlsx", sheet_name="BIOLOG"
)

df_cleaner = Cleaner()
values_df = df_cleaner.cleaNorganize(biolog_values, compounds=biolog_compounds)
biolog_raw_id_df = df_cleaner.cleaNorganize(biolog_ids, compounds=biolog_compounds)
df_535 = df_cleaner.cleaNorganize(abs_535, compounds=biolog_compounds)
df_575 = df_cleaner.cleaNorganize(abs_575, compounds=biolog_compounds)
df_590 = df_cleaner.cleaNorganize(abs_590, compounds=biolog_compounds)
df_750 = df_cleaner.cleaNorganize(abs_750, compounds=biolog_compounds)

########################

# * BEAT BWA
beat_values = BEAT(df=values_df)
beat_values.calculation_dif()

# * BEAT 590
beat_590 = BEAT(df=df_590)
beat_590.calculation_dif()

########################

bnt_values = BNT(dataframe=values_df)
bnt_535 = BNT(dataframe=df_535)
bnt_575 = BNT(dataframe=df_575)
bnt_590 = BNT(dataframe=df_590)
bnt_750 = BNT(dataframe=df_750)

##########################
# * Strain IDs Dataframe

strain_ids = reader.find_dataframe("strain_ids.xlsx", sheet_name="IDs")
strain_ids.set_index("ID-number", inplace=True)

##########################

biolog_ids_test = reader.find_dataframe(
    "plate_readout_categorizations_test.xlsx",
    sheet_name="BIOLOG",
)
biolog_ids_test = df_cleaner.cleaNorganize(biolog_ids_test, compounds=biolog_compounds)

##########################

# * BNT BWA
bnt_values.identify_result_type(
    dataframe=bnt_values.complete_df,
    limit_i=1.0,
    limit_b=1.5,
)

# * BNT 590
bnt_590.identify_result_type(
    dataframe=bnt_590.normalized_df,
    limit_i=1.0,
    limit_b=1.5,
)

##########################

baa_test = BAA()

# * Pairgrids for BWA data, both raw and normalized

test = baa_test.distribution_analysis(
    id_df=biolog_ids_test,
    tgt_df=beat_values.dif_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=False,
)

test.sort_values("Range Dif").reset_index(drop=True)

baa_test.distribution_analysis(
    id_df=bnt_values.id_df,
    tgt_df=beat_values.dif_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=False,
)

baa_test.distribution_analysis(
    id_df=biolog_ids_test,
    tgt_df=bnt_values.normalized_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=False,
)

baa_test.distribution_analysis(
    id_df=bnt_values.id_df,
    tgt_df=bnt_values.normalized_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=False,
)

##########################

# * Pairgrids for 590 data, both raw and normalized

abs_test = baa_test.distribution_analysis(
    id_df=biolog_ids_test,
    tgt_df=beat_590.dif_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=True,
)

abs_test.sort_values("Range Dif").reset_index(drop=True)


baa_test.distribution_analysis(
    id_df=bnt_590.id_df,
    tgt_df=beat_590.dif_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=True,
)

abs_test2 = baa_test.distribution_analysis(
    id_df=biolog_ids_test,
    tgt_df=bnt_590.normalized_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=True,
)

abs_test2.sort_values("Range Dif").reset_index(drop=True)

baa_test.distribution_analysis(
    id_df=bnt_590.id_df,
    tgt_df=bnt_590.normalized_df,
    strain_id_df=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    diag_cut=2,
    save_figs=True,
)
