import set_cwd

from notebooks import biolog_variables as bv
from notebooks import exploratory_variables as ev

from src import Reader, Cleaner, BEAT, BNT, BAA

##########################

# import pandas_settings
from src import PathOrganizer
from src import Model_Exploration_Tool as es
from src import ExploratoryAid as EA

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import cobra
import re

# %matplotlib inline

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)

po = PathOrganizer()

baa_test = BAA()
test = EA()

%%time
biomass_mAssay_results, ATPM_mAssay_results = test.metabolite_assay(
    target_models=ev.draftModels,
    medium=ev.kwoji_updated,
    closed=ev.closed_uptake,
    essential=ev.essential,
    metabolites=ev.biolog_met_df,
    medium_df=ev.kwoji_met_df,
    optimizationType="llFBA",
)

test.save_dataframe_to_pickle(dataframe=biomass_mAssay_results)
test.save_dataframe_to_pickle(dataframe=ATPM_mAssay_results)

biomass_mAssay_results = test.load_dataframe_from_pickle()
ATPM_mAssay_results = test.load_dataframe_from_pickle()

##########################

baa_test.generate_SummaryCategoriesMetAssayTable(
    categoriesDataframe_BIOLOG=bv.biolog_ids_test,
    categoriesDataframe_BWA=bv.bnt_values.id_df,
    categoriesDataframe_590nm=bv.bnt_590.id_df,
    categoriesDataframe_750nm=bv.bnt_750.id_df,
    strainInfoDataframe=bv.strain_ids,
    biomass_metAssayData=biomass_mAssay_results,
    mainAPC_metAssayData=ATPM_mAssay_results,
    saveFigures=False,
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

## BWA, BIOLOG
baa_test.generate_summaryPieChartReports(
    targetData=bv.beat_values.dif_df,
    categoriesDataframe=bv.biolog_ids_test,
    strainIDsDataframe=bv.strain_ids,
    metabolitesDataframe=bv.mets_dataframe,
    saveFigures=False,
)

## BWA, BNT
baa_test.generate_summaryPieChartReports(
    targetData=bv.beat_values.dif_df,
    categoriesDataframe=bv.bnt_values.id_df,
    strainIDsDataframe=bv.strain_ids,
    metabolitesDataframe=bv.mets_dataframe,
    saveFigures=False,
)

## 590nm, BIOLOG
baa_test.generate_summaryPieChartReports(
    targetData=bv.beat_590.dif_df,
    categoriesDataframe=bv.biolog_ids_test,
    strainIDsDataframe=bv.strain_ids,
    metabolitesDataframe=bv.mets_dataframe,
    saveFigures=False,
)

## 590nm, BNT
baa_test.generate_summaryPieChartReports(
    targetData=bv.beat_590.dif_df,
    categoriesDataframe=bv.bnt_590.id_df,
    strainIDsDataframe=bv.strain_ids,
    metabolitesDataframe=bv.mets_dataframe,
    saveFigures=False,
)

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
