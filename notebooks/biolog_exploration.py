import set_cwd
import pandas as pd
import numpy as np

pd.set_option("display.precision", 3)
pd.set_option("display.max_columns", 200)
pd.set_option("display.max_rows", 200)

from src import Reader, Cleaner, BEAT, BNT

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

bnt_values = BNT(dataframe=values_df)
bnt_535 = BNT(dataframe=df_535)
bnt_575 = BNT(dataframe=df_575)
bnt_590 = BNT(dataframe=df_590)
bnt_750 = BNT(dataframe=df_750)

# ----------------------------------------------------------------------------------- #

# beat_values = BEAT(df=values_df)
# beat_values.data_analysis()
# beat_values.display_categories(dataframe=beat_values.complete_df)
# res_values = beat_values.return_categorized_data(dataframe=beat_values.complete_df)

# beat_values.stats_df["23-47"]
# beat_values.complete_df["23-42"]

# beat_535 = BEAT(df=df_535)
# beat_535.data_analysis()
# beat_535.display_categories(dataframe=beat_535.complete_df)
# res_535 = beat_535.return_categorized_data(dataframe=beat_535.complete_df)

# beat_575 = BEAT(df=df_575)
# beat_575.data_analysis()
# beat_575.display_categories(dataframe=beat_575.complete_df)
# res_575 = beat_575.return_categorized_data(dataframe=beat_575.complete_df)

# beat_590 = BEAT(df=df_590)
# beat_590.data_analysis()
# beat_590.display_categories(dataframe=beat_590.complete_df)
# res_590 = beat_590.return_categorized_data(dataframe=beat_590.complete_df)

# beat_590.stats_df["23-42"]

# beat_590.stats_df["23-47"]

# beat_750 = BEAT(df=df_750)
# beat_750.data_analysis()
# beat_750.display_categories(dataframe=beat_750.complete_df)
# res_750 = beat_750.return_categorized_data(dataframe=beat_750.complete_df)

# res_values.to_excel("res_values.xlsx", engine="openpyxl")
# res_535.to_excel("res_535.xlsx", engine="openpyxl")
# res_575.to_excel("res_575.xlsx", engine="openpyxl")
# res_590.to_excel("res_590.xlsx", engine="openpyxl")
# res_750.to_excel("res_750.xlsx", engine="openpyxl")

# ----------------------------------------------------------------------------------- #

# bnt_values.yeojohn_normalization(dataframe=bnt_values.complete_df)

# bnt_values.generate_skewness_data()

# bnt_values.skew_df.loc[["23-36", "23-38", "23-47"]].map("{:.2f}".format)

# lala = bnt_values.identify_result_type(dataframe=bnt_values.normalized_df, limit_i=1.5, limit_b=2.0)
# lala

# bnt_values.yeojohn_normalization(dataframe=bnt_values.complete_df)
# coco = bnt_values.normalized_df["23-47"]
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

# biolog_raw_id_df

# normalized_df = bnt_values.yeojohn_normalization(dataframe=bnt_values.complete_df["23-47"])

# test_excel = bnt_values.display_categories(dataframe=normalized_df, limit_i=1.5, limit_b=2.0)

# ad_test = BEAT(df=values_df)

# ad_test.calculation_dif()

# ad_test.data_analysis()

# ad_test.display_categories(
#     dataframe=ad_test.complete_df, tag_df=ad_test.id_df, columns="23-47"
# )

# lolo = ad_test.return_categorized_data(dataframe=ad_test.complete_df, columns="23-47")

# lolo.to_excel("testexcel.xlsx", engine="openpyxl")

# ----------------------------------------------------------------------------------- #

# cat_val = bnt_values.display_categories(
#     dataframe=bnt_values.complete_df["23-47"], limit_i=1.5, limit_b=2.0
# )
# cat_535 = bnt_535.display_categories(
#     dataframe=bnt_535.complete_df["23-42"], limit_i=1.5, limit_b=2.0
# )
# cat_575 = bnt_575.display_categories(
#     dataframe=bnt_575.complete_df["23-42"], limit_i=1.5, limit_b=2.0
# )
# cat_590 = bnt_590.display_categories(
#     dataframe=bnt_590.complete_df["23-42"], limit_i=1.5, limit_b=2.0
# )
# cat_750 = bnt_750.display_categories(
#     dataframe=bnt_750.complete_df["23-42"], limit_i=1.5, limit_b=2.0
# )

# cat_val.to_excel("cat_val.xlsx", engine="openpyxl")
# cat_535.to_excel("cat_535.xlsx", engine="openpyxl")
# cat_575.to_excel("cat_575.xlsx", engine="openpyxl")
# cat_590.to_excel("cat_590.xlsx", engine="openpyxl")
# cat_750.to_excel("cat_750.xlsx", engine="openpyxl")

# bnt_values.display_yeojohn_norm(target_columns="23-42")
# bnt_values.display_yeojohn_norm(target_columns="23-47")
# bnt_535.display_yeojohn_norm(target_columns="23-42")
# bnt_575.display_yeojohn_norm(target_columns="23-42")
# bnt_590.display_yeojohn_norm(target_columns="23-47")
# bnt_750.display_yeojohn_norm(target_columns="23-42")

# bnt_values.generate_skewness_data()
# bnt_535.generate_skewness_data()
# bnt_575.generate_skewness_data()
# bnt_590.generate_skewness_data()
# bnt_750.generate_skewness_data()

#####################################

test_df = bnt_values.identify_best_boundaries(
    cat_df=biolog_ids, target_df=bnt_values.complete_df
)

test_df

###############################

import set_cwd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib as Path

pd.set_option("display.precision", 0)

from src import Reader, Cleaner, BEAT, BNT

reader = Reader()

maja_df = reader.find_dataframe(
    "output_df_all_strains.xlsx", sheet_name="Sheet1"
).set_index("Compounds")

strains_list = [50, 49, 48, 47, 46, 45, 44, 43, 42, 40, 39, 38, 37, 36, 34, 33, 32]

####

from datetime import datetime
from pathlib import Path


# Generate a new folder named by Year_Month_Day if not available, get its path
def generate_daily_folder() -> Path:

    # Establish today's date
    today = datetime.today()

    # Define as string
    datestr = today.strftime("%Y_%m_%d")

    # Define and resolve the path to the new data folder
    date_path = Path.cwd() / "reports" / "data" / datestr
    date_path.resolve()

    # Generate folder
    Path(date_path).mkdir(parents=True, exist_ok=True)

    # Return folder's path as a Path object
    return date_path


####

today_data_fdr = generate_daily_folder()

############################

for col in maja_df.columns:

    maja_df[col] = maja_df[col] - maja_df.loc["Water", col]

for sample_num in strains_list:

    all_data_one_strain = [col for col in maja_df if col.startswith(str(sample_num))]

    cat_all_plates = [col for col in all_data_one_strain if "Cat" in col]
    bwa_all_plates = [col for col in all_data_one_strain if "BWA" in col]
    abs_all_plates = [col for col in all_data_one_strain if "590_nm" in col]
    growth_all_plates = [col for col in all_data_one_strain if "750_nm" in col]

    fig, axes = plt.subplots(1, 3, figsize=(30, 6))

    fig.suptitle(f"BWA for Sample {sample_num}", fontsize=20)

    plt.subplot(1, 3, 1)
    sns.histplot(
        ax=axes[0],
        data=maja_df,
        x=bwa_all_plates[0],
        kde=False,
        hue=cat_all_plates[0],
        multiple="stack",
        bins=30,
        legend=False,
    )
    plt.title("Plate 1", fontsize=15)
    plt.xlabel("BWA Value")

    plt.subplot(1, 3, 2)
    sns.histplot(
        ax=axes[1],
        data=maja_df,
        x=bwa_all_plates[1],
        kde=False,
        hue=cat_all_plates[1],
        multiple="stack",
        bins=30,
        legend=False,
    )
    plt.title("Plate 2", fontsize=15)
    plt.xlabel("BWA Value")

    plt.subplot(1, 3, 3)
    sns.histplot(
        ax=axes[2],
        data=maja_df,
        x=bwa_all_plates[2],
        kde=False,
        hue=cat_all_plates[2],
        multiple="stack",
        bins=30,
        legend=False,
    )
    plt.title("Plate 3", fontsize=15)
    plt.xlabel("BWA Value")

    img_path = today_data_fdr / f"BWA_{sample_num}"
    img_path.resolve()

    plt.savefig(img_path)
    plt.show

##################################

strains = [47]

for strain in strains:

    all_data_one_strain = [col for col in maja_df if col.startswith(str(strain))]

    # Separate all the different data types into specific lists
    bwa_all_plates = [col for col in all_data_one_strain if "BWA" in col]
    abs_all_plates = [col for col in all_data_one_strain if "nm" in col]
    biolog_all_plates = [col for col in all_data_one_strain if "590_nm" in col]
    growth_all_plates = [col for col in all_data_one_strain if "750_nm" in col]

    cat_all_plates = [col for col in all_data_one_strain if "Cat" in col]

    # Calculate the mean values for the categorization values
    means = (
        maja_df.loc[:, cat_all_plates[0]]
        + maja_df.loc[:, cat_all_plates[1]]
        + maja_df.loc[:, cat_all_plates[2]]
    ) / 3

    # Assign a general mean value based on the calculated mean
    cat_mean = []
    for mean in means:
        if mean < 0.5:
            cat_mean.append(0)
        elif mean < 1.5:
            cat_mean.append(1)
        else:
            cat_mean.append(2)

    mean_df = maja_df.copy()
    mean_df["Cat_mean_" + str(strain)] = cat_mean

    dataf_ab = maja_df.loc[:, abs_all_plates]

    cols = dataf_ab.columns.tolist()
    cols = cols[::2] + cols[1::2]
    dataf_ab = dataf_ab[cols]

    # Generate the plots
    ## Structure of the plots
    l_split = [[0, 48], [48, 96]]

    ## Limit value set for the for the range for the y axis
    ylimit = dataf_ab.values.max() + 0.05

    for n in range(len(l_split)):

        split = l_split[n]

        ax = dataf_ab.iloc[split[0] : split[1], :].plot(
            kind="bar",
            legend=False,
            width=0.7,
            color=["#33adff", "#3333ff", "#00004d", "#ffcc00", "#cc6600", "#991f00"],
            figsize=(25, 4),
        )

        ax.set_ylim(0, ylimit)

        for i, value in enumerate(cat_mean[split[0] : split[1]]):
            if value == 2:
                # This are the compounds considered by BIOLOG to be positive
                ax.get_xticklabels()[i].set_color("green")
            if value == 1:
                # This are the compounds considered by BIOLOG to be boundary
                ax.get_xticklabels()[i].set_color("orange")

        plt.xticks(rotation=45, ha="right")
        # plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)

        plt.show()

##################################

# strains = [50, 49, 48, 47, 46, 45, 44, 43, 42, 40, 39, 38, 37, 36, 34, 33, 32]

import re

new_strains = bnt_values.complete_df.columns.get_level_values(level=0).unique().tolist()

pd.set_option("display.precision", 3)

strains = []
for col in new_strains:

    x = re.sub("^23-", "", col)
    strains.append(x)

for strain in strains:

    all_data_one_strain = [col for col in maja_df if col.startswith(str(strain))]

    # Separate all the different data types into specific lists
    bwa_all_plates = [col for col in all_data_one_strain if "BWA" in col]
    abs_all_plates = [col for col in all_data_one_strain if "nm" in col]
    biolog_all_plates = [col for col in all_data_one_strain if "590_nm" in col]
    growth_all_plates = [col for col in all_data_one_strain if "750_nm" in col]

    cat_all_plates = [col for col in all_data_one_strain if "Cat" in col]

    # Calculate the mean values for the categorization values
    means = (
        maja_df.loc[:, cat_all_plates[0]]
        + maja_df.loc[:, cat_all_plates[1]]
        + maja_df.loc[:, cat_all_plates[2]]
    ) / 3

    # Assign a general mean value based on the calculated mean
    cat_mean = []
    for mean in means:
        if mean < 0.5:
            cat_mean.append(0)
        elif mean < 1.5:
            cat_mean.append(1)
        else:
            cat_mean.append(2)

    mean_df = maja_df.copy()
    mean_df["Cat_mean_" + str(strain)] = cat_mean

    mean_all_plates = [col for col in mean_df if f"Cat_mean_{str(strain)}" in col]

    target_all_plates = biolog_all_plates + mean_all_plates

    target_df = mean_df.loc[:, target_all_plates]

    g = sns.PairGrid(
        target_df,
        hue=f"Cat_mean_{str(strain)}",
        corner=False,
        palette="colorblind",
        height=4,
    )

    g.map_diag(sns.histplot, multiple="stack", element="step")
    g.map_upper(sns.scatterplot)
    g.map_lower(sns.kdeplot, levels=10, thresh=0.05)

    for ax in g.axes.flat:

        xy_min = mean_df.loc[:, biolog_all_plates].min().min()

        xy_max = mean_df.loc[:, biolog_all_plates].max().max()

        xy_min = xy_min - (0.25 * xy_max)
        xy_max = xy_max + (0.25 * xy_max)

        ax.set_xlim([xy_min, xy_max])
        ax.set_ylim([xy_min, xy_max])

    g.add_legend()

    img_path = today_data_fdr / f"Strain_{str(strain)}_DistComp"
    img_path.resolve()

    g.savefig(img_path)

    plt.show()

######################
