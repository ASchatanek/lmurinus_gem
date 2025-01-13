######################

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

###############################
