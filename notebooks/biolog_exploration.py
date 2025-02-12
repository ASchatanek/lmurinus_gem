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
biolog_values = reader.find_dataframe(
    file_name="plate_readout_numbers.xlsx",
    sheet_name="BIOLOG",
)
abs_535 = reader.find_dataframe(
    file_name="plate_readout_numbers.xlsx",
    sheet_name="Abs_535",
)
abs_575 = reader.find_dataframe(
    file_name="plate_readout_numbers.xlsx",
    sheet_name="Abs_575",
)
abs_590 = reader.find_dataframe(
    file_name="plate_readout_numbers.xlsx",
    sheet_name="Abs_590",
)
abs_750 = reader.find_dataframe(
    file_name="plate_readout_numbers.xlsx",
    sheet_name="Abs_750",
)

biolog_compounds_df = reader.find_dataframe(
    file_name="plate_readout_numbers.xlsx",
    sheet_name="Compounds",
)

biolog_compounds = reader.get_compounds_list(biolog_compounds_df)

biolog_ids = reader.find_dataframe(
    file_name="plate_readout_categorizations.xlsx",
    sheet_name="BIOLOG",
)

########################
# * Clean and Organize

df_cleaner = Cleaner()
values_df = df_cleaner.cleaNorganize(
    dataframe=biolog_values,
    compounds=biolog_compounds,
)
biolog_raw_id_df = df_cleaner.cleaNorganize(
    dataframe=biolog_ids,
    compounds=biolog_compounds,
)
df_535 = df_cleaner.cleaNorganize(
    dataframe=abs_535,
    compounds=biolog_compounds,
)
df_575 = df_cleaner.cleaNorganize(
    dataframe=abs_575,
    compounds=biolog_compounds,
)
df_590 = df_cleaner.cleaNorganize(
    dataframe=abs_590,
    compounds=biolog_compounds,
)
df_750 = df_cleaner.cleaNorganize(
    dataframe=abs_750,
    compounds=biolog_compounds,
)

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

strain_ids = reader.find_dataframe(
    file_name="strain_ids.xlsx",
    sheet_name="IDs",
)
strain_ids.set_index("ID-number", inplace=True)

##########################

biolog_ids_test = reader.find_dataframe(
    file_name="plate_readout_categorizations_test.xlsx",
    sheet_name="BIOLOG",
)
biolog_ids_test = df_cleaner.cleaNorganize(
    dataframe=biolog_ids_test,
    compounds=biolog_compounds,
)

mets_dataframe = reader.find_dataframe(
    file_name="plate_readout_categorizations_test.xlsx",
    sheet_name="Compounds",
)

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

## BWA, BIOLOG
baa_test.generate_summary_report(
    data_dataframe=beat_values.dif_df,
    categories_dataframe=biolog_ids_test,
    strain_id_dataframe=strain_ids,
    mets_dataframe=mets_dataframe,
    save_figures=True,
)

## BWA, BNT
baa_test.generate_summary_report(
    data_dataframe=beat_values.dif_df,
    categories_dataframe=bnt_values.id_df,
    strain_id_dataframe=strain_ids,
    mets_dataframe=mets_dataframe,
    save_figures=True,
)

## 590nm, BIOLOG
baa_test.generate_summary_report(
    data_dataframe=beat_590.dif_df,
    categories_dataframe=biolog_ids_test,
    strain_id_dataframe=strain_ids,
    mets_dataframe=mets_dataframe,
    save_figures=True,
)

##590nm, BNT
baa_test.generate_summary_report(
    data_dataframe=beat_590.dif_df,
    categories_dataframe=bnt_590.id_df,
    strain_id_dataframe=strain_ids,
    mets_dataframe=mets_dataframe,
    save_figures=True,
)

# * Pairgrids for BWA data, both raw and normalized

## BWA (raw), BIOLOG
baa_test.values_distribution(
    data_dataframe=beat_values.dif_df,
    categories_dataframe=biolog_ids_test,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

## BWA (raw), BNT
baa_test.values_distribution(
    data_dataframe=beat_values.dif_df,
    categories_dataframe=bnt_values.id_df,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

## BWA (normalized), BIOLOG
baa_test.values_distribution(
    data_dataframe=bnt_values.normalized_df,
    categories_dataframe=biolog_ids_test,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

## BWA (normalized), BNT
baa_test.values_distribution(
    data_dataframe=bnt_values.normalized_df,
    categories_dataframe=bnt_values.id_df,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

##########################

# * Pairgrids for 590 data, both raw and normalized

## 590nm (raw), BIOLOG
baa_test.values_distribution(
    data_dataframe=beat_590.dif_df,
    categories_dataframe=biolog_ids_test,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)

## 590nm (raw), BNT
baa_test.values_distribution(
    data_dataframe=beat_590.dif_df,
    categories_dataframe=bnt_590.id_df,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)

## 590nm (normalized), BIOLOG
abs_test2 = baa_test.values_distribution(
    data_dataframe=bnt_590.normalized_df,
    categories_dataframe=biolog_ids_test,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)

## 590nm (normalized), BNT
baa_test.values_distribution(
    data_dataframe=bnt_590.normalized_df,
    categories_dataframe=bnt_590.id_df,
    strain_id_dataframe=strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)
