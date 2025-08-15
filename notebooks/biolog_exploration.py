import re
import cobra
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import seaborn as sns

from notebooks import biolog_variables as bv
from notebooks import exploratory_variables as ev

##########################
# import pandas_settings
from src import BAA, BEAT, BNT, Cleaner, PathOrganizer, Reader
from src import ExploratoryAid as EA
from src import Model_Exploration_Tool as es

# %matplotlib inline

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)

po = PathOrganizer()

baa_test = BAA()
test = EA()

okabe_ito_palette = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#000000",  # black
    "#BBBBBB",  # light grey
    "#7E2954",  # purple
]


##########################

summary = baa_test.testing_something(
    categoriesDataframe_BIOLOG=bv.biolog_ids_test,
    categoriesDataframe_BWA=bv.bnt_values.id_df,
    categoriesDataframe_590nm=bv.bnt_590.id_df,
)

summary = summary.reorder_levels(["Sample", "TP", "Data Types"], axis=1).sort_index(
    level=0, axis=1
)

result = (
    summary.T.groupby(
        ["Sample", "TP"],
        group_keys=True,
    )
    .mean()
    .round(decimals=2)
)

okay = result.sum(axis=1)
aver = okay.loc[okay.groupby(["Sample"]).idxmax()]

aver = list(aver.index)

bv.strain_ids

target_results = result.loc[aver, :]

strains = result.index.get_level_values(level=0).unique().tolist()
strain_names = bv.strain_ids.loc[strains, "Strain"].tolist()

###

okabe_ito_rgb = [
    colors.to_rgb("#DDDDDD"),  # light grey
    colors.to_rgb("#7E2954"),  # purple
    colors.to_rgb("#009E73"),  # bluish green
]

# Create arrays repeating each color to match original segment lengths
colormap_under_10 = np.tile(okabe_ito_rgb[0], (49, 1))  # light grey
colormap_10_30 = np.tile(okabe_ito_rgb[1], (33, 1))  # sky blue
colormap_above_30 = np.tile(okabe_ito_rgb[2], (18, 1))  # bluish green

# Combine all into one array
colors_combined = np.vstack((colormap_under_10, colormap_10_30, colormap_above_30))
mymap = colors.LinearSegmentedColormap.from_list("colormap_combined", colors_combined)

plt.figure(figsize=(20, 4), dpi=600)

sns.set_theme(font_scale=0.75)
sns.heatmap(
    data=target_results,
    vmin=0,
    vmax=2,
    cmap=mymap,
    linewidths=0.75,
    yticklabels=strain_names,
)
plt.ylabel("Strains")
plt.xlabel("Metabolites")
# plt.xticks(rotation=45, ha="right")

plt.show()


###################

met_groups = ev.biolog_met_complete_df.loc[:, ("Metabolite", "Group")]
met_groups = met_groups.set_index("Metabolite")


target_results = result.loc[aver, :]

targetresults_newIndex = target_results.index.droplevel(level=1)
target_results.index = targetresults_newIndex

group_list = met_groups["Group"].to_list()

target_results.loc["Group"] = group_list

martin = target_results.loc[target_results.index != "Group", :].astype(np.float64)
martin.dtypes
martin = martin.astype(np.float64)

martin.dtypes
target_results.dtypes

SPECIES_COLORS = target_results.loc["Group", :].map(
    {
        "Miscellaneous / Other": "#000000",
        "Carbohydrates & Derivatives": "#D55E00",
        "Organic Acids & Derivatives": "#56B4E9",
        "Amino Acids & Peptides": "#CC79A7",
        "Nucleosides & Nucleotides": "#0072B2",
    }
)

met_groups = ev.biolog_met_complete_df.loc[:, ("Metabolite", "Group")]
met_groups = met_groups.set_index("Metabolite")

group_list = met_groups["Group"].to_list()

Group_colors = {
    "Miscellaneous / Other": "#000000",
    "Carbohydrates & Derivatives": "#D55E00",
    "Organic Acids & Derivatives": "#56B4E9",
    "Amino Acids & Peptides": "#CC79A7",
    "Nucleosides & Nucleotides": "#0072B2",
}

metGroup_colors = [Group_colors[x] for x in group_list]

###################

sns.set_theme(rc={"figure.dpi": 600})
sns.set_theme(rc={"xtick.bottom": True, "ytick.left": True})
sns.set_theme(font_scale=0.75)
g = sns.clustermap(
    data=martin,
    figsize=(17, 6),
    annot=False,
    dendrogram_ratio=(0.05, 0.2),
    col_cluster=False,
    yticklabels=strain_names,
    row_cluster=False,
    cmap=mymap,
    linewidth=0.75,
    col_colors=metGroup_colors,
)

g.ax_cbar.remove()
g.ax_heatmap.yaxis.set_ticks_position("left")
g.ax_heatmap.yaxis.set_label_position("left")

col_pos = g.ax_col_colors.get_position()
heatmap_pos = g.ax_heatmap.get_position()

# Shift col_colors axis below heatmap
g.ax_col_colors.set_position(
    [
        heatmap_pos.x0,
        heatmap_pos.y0 - 0.01,  # small offset
        heatmap_pos.width,
        col_pos.height,
    ]
)

# Adjust heatmap to shift it up if necessary
g.ax_heatmap.set_position(
    [
        heatmap_pos.x0,
        heatmap_pos.y0 + 0.01,
        heatmap_pos.width,
        heatmap_pos.height,
    ]
)
# g.ax_heatmap.xaxis.set_tick_params(which="minor")

g.ax_heatmap.set_xlabel(
    "Metabolites",
    fontdict={
        "size": 12,
        "weight": "bold",
    },
)
g.ax_heatmap.set_ylabel(
    "Strains",
    fontdict={
        "size": 12,
        "weight": "bold",
    },
)

plt.show()

##########################

crucial, variable = test.generate_essentialMetabolites_dictionaries(
    targetModels=ev.draftModels,
    mediumMetabolites=ev.kwoji_medium,
    mediumMetabolitesDataframe=ev.kwoji_medium_df,
    closedMetabolites=ev.closed_uptake,
)

best_eMets = test.essentialMetabolites_report(
    targetModels=ev.draftModels,
    mediumMetabolites=ev.kwoji_medium,
    mediumMetabolitesDataframe=ev.kwoji_medium_df,
    closedMetabolites=ev.closed_uptake,
    crucial_eMetabolites=crucial,
    variable_eMetabolites=variable,
)

biomass_mAssay_results, ATPM_mAssay_results, metabolitesExistence = (
    test.metabolite_assay(
        targetModels=ev.draftModels,
        mediumMetabolites=ev.kwoji_medium,
        closedMetabolites=ev.closed_uptake,
        targetModels_EssentialMetabolites=best_eMets,
        assayMetabolitesDataframe=ev.biolog_met_df,
        mediumMetabolitesDataframe=ev.kwoji_medium_df,
        optimizationType="FBA",
    )
)

pd.concat(
    [biomass_mAssay_results, ATPM_mAssay_results, metabolitesExistence],
    axis=1,
    keys=["bio", "ATPM", "met"],
).sort_index(axis=1, level=1)

##########################

baa_test.generate_SummaryCategoriesMetAssayTable(
    categoriesDataframe_BIOLOG=bv.biolog_ids_test,
    categoriesDataframe_BWA=bv.bnt_values.id_df,
    categoriesDataframe_590nm=bv.bnt_590.id_df,
    categoriesDataframe_750nm=bv.bnt_750.id_df,
    strainInfoDataframe=bv.strain_ids,
    biomass_metAssayData=biomass_mAssay_results,
    mainATP_metAssayData=ATPM_mAssay_results,
    saveFigures=True,
)

baa_test.generate_SummaryCategoriesMetAssayTable(
    categoriesDataframe_BIOLOG=bv.biolog_ids_test,
    categoriesDataframe_BWA=bv.bnt_values.id_df,
    categoriesDataframe_590nm=bv.bnt_590.id_df,
    strainInfoDataframe=bv.strain_ids,
    Biomass_metAssayData=ATPM_mAssay_results,
    saveFigures=False,
)

##########################

# * Pairgrids for BWA data, both raw and normalized

## BWA (raw), BIOLOG
baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.beat_values.dif_df,
    categoriesDataframe=bv.biolog_ids_test,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)

## BWA (raw), BNT
baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.beat_values.dif_df,
    categoriesDataframe=bv.bnt_values.id_df,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)

## BWA (normalized), BIOLOG
baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.bnt_values.normalized_df,
    categoriesDataframe=bv.biolog_ids_test,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)

## BWA (normalized), BNT
baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.bnt_values.normalized_df,
    categoriesDataframe=bv.bnt_values.id_df,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)

##########################

# * Pairgrids for 590 data, both raw and normalized

## 590nm (raw), BIOLOG
baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.beat_590.dif_df,
    categoriesDataframe=bv.biolog_ids_test,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)

## 590nm (raw), BNT
baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.beat_590.dif_df,
    categoriesDataframe=bv.bnt_590.id_df,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)

## 590nm (normalized), BIOLOG
abs_test2 = baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.bnt_590.normalized_df,
    categoriesDataframe=bv.biolog_ids_test,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)

## 590nm (normalized), BNT
baa_test.generate_plateDistributionsPairGrids(
    targetData=bv.bnt_590.normalized_df,
    categoriesDataframe=bv.bnt_590.id_df,
    strainIDsDataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    saveFigures=False,
)
