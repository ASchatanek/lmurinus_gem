import set_cwd
import pandas_settings
import pandas as pd
import numpy as np

from src import Reader, Cleaner, BAT, BNT

reader = Reader()

biolog_values = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="BIOLOG")
abs_values = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="Abs_535")
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
abs_df = df_cleaner.cleaNorganize(abs_values, compounds=biolog_compounds)

cc = BAT(values_df)
cc.data_analysis()

cc2 = BAT(abs_df)
cc2.data_analysis()

test = BNT(dataframe=values_df)

test2 = BNT(dataframe=abs_df)

df = test.yeojohn_normalization().copy()

std = df.describe().copy()

std = std.loc["std", :]

test.display_data_categorization(dataframe=df["23-47"])

ef = test2.yeojohn_normalization().copy()

test2.display_data_categorization(dataframe=ef)
