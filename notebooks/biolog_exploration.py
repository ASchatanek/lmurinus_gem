import biolog_variables as bv

from src import Reader, Cleaner, BEAT, BNT, BAA

##########################

baa_test = BAA()

baa_test.generate_super_report(
    cat_biolog_df=bv.biolog_ids_test,
    cat_bwa_norm_df=bv.bnt_values.id_df,
    cat_590nm_norm_df=bv.bnt_590.id_df,
    strain_id_df=bv.strain_ids,
)

## BWA, BIOLOG
baa_test.generate_summary_report(
    data_dataframe=bv.beat_values.dif_df,
    categories_dataframe=bv.biolog_ids_test,
    strain_id_dataframe=bv.strain_ids,
    mets_dataframe=bv.mets_dataframe,
    save_figures=False,
)

## BWA, BNT
baa_test.generate_summary_report(
    data_dataframe=bv.beat_values.dif_df,
    categories_dataframe=bv.bnt_values.id_df,
    strain_id_dataframe=bv.strain_ids,
    mets_dataframe=bv.mets_dataframe,
    save_figures=False,
)

## 590nm, BIOLOG
baa_test.generate_summary_report(
    data_dataframe=bv.beat_590.dif_df,
    categories_dataframe=bv.biolog_ids_test,
    strain_id_dataframe=bv.strain_ids,
    mets_dataframe=bv.mets_dataframe,
    save_figures=False,
)

## 590nm, BNT
baa_test.generate_summary_report(
    data_dataframe=bv.beat_590.dif_df,
    categories_dataframe=bv.bnt_590.id_df,
    strain_id_dataframe=bv.strain_ids,
    mets_dataframe=bv.mets_dataframe,
    save_figures=False,
)

# * Pairgrids for BWA data, both raw and normalized

## BWA (raw), BIOLOG
baa_test.values_distribution(
    data_dataframe=bv.beat_values.dif_df,
    categories_dataframe=bv.biolog_ids_test,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

## BWA (raw), BNT
baa_test.values_distribution(
    data_dataframe=bv.beat_values.dif_df,
    categories_dataframe=bv.bnt_values.id_df,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

## BWA (normalized), BIOLOG
baa_test.values_distribution(
    data_dataframe=bv.bnt_values.normalized_df,
    categories_dataframe=bv.biolog_ids_test,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

## BWA (normalized), BNT
baa_test.values_distribution(
    data_dataframe=bv.bnt_values.normalized_df,
    categories_dataframe=bv.bnt_values.id_df,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=False,
)

##########################

# * Pairgrids for 590 data, both raw and normalized

## 590nm (raw), BIOLOG
baa_test.values_distribution(
    data_dataframe=bv.beat_590.dif_df,
    categories_dataframe=bv.biolog_ids_test,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)

## 590nm (raw), BNT
baa_test.values_distribution(
    data_dataframe=bv.beat_590.dif_df,
    categories_dataframe=bv.bnt_590.id_df,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)

## 590nm (normalized), BIOLOG
abs_test2 = baa_test.values_distribution(
    data_dataframe=bv.bnt_590.normalized_df,
    categories_dataframe=bv.biolog_ids_test,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)

## 590nm (normalized), BNT
baa_test.values_distribution(
    data_dataframe=bv.bnt_590.normalized_df,
    categories_dataframe=bv.bnt_590.id_df,
    strain_id_dataframe=bv.strain_ids,
    kde_levels=5,
    kde_thresh=0.1,
    save_figures=True,
)
