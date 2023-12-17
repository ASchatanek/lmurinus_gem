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

test = BNT(dataframe=values_df)
test2 = BNT(dataframe=abs_df)

test.yeojohn_normalization(dataframe=test.complete_df)

test.generate_skewness_data(dataframe=test.complete_df)

test.skew_df.loc[["23-36", "23-38", "23-47"]].map("{:.2f}".format)

lala = test.identify_result_type(dataframe=test.normalized_df, limit_i=1.5, limit_b=2.0)
lala

# test.yeojohn_normalization(dataframe=test.complete_df)
# coco = test.normalized_df["23-47"]
# coco.describe().map("{:.2f}".format)
# coco = coco.rename_axis(["Timepoints", "Plate Number"], axis="columns")

# import matplotlib.pyplot as plt
# import seaborn as sns

# columns = coco.columns.get_level_values(level=0).unique()
# for column in columns:
#     plt.figure(figsize=(8, 8))  # Set the figure size (optional)
#     sns.boxplot(data=coco[column], width=0.5, whis=1)
#     plt.title(f"Boxplots of Plates at Timepoint {column}")  # Set title with column name
#     plt.xlabel("Plates")  # Set x-axis label (if needed)
#     plt.ylabel("Normalized Data")  # Set y-axis label (if needed)
#     plt.grid(True)  # Add gridlines (optional)
#     plt.show()

biolog_raw_id_df
