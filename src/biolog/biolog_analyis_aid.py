import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy import stats

from IPython.display import display


class BAA:

    def __init__(self):
        pass

    def iterate(self, tgt_df: pd.DataFrame):
        samples = tgt_df.columns.get_level_values(level=0).unique()

        for sample in samples:

            tps = tgt_df.loc[:, sample].columns.get_level_values(level=0).unique()

            for tp in tps:

                plates = (
                    tgt_df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )

                for plate in plates:
                    tgt = (sample, tp, plate)
                    yield (tgt)

    def transform_ids_to_num(
        self,
        categories_dataframe: pd.DataFrame,
    ) -> pd.DataFrame:

        cat_df = categories_dataframe

        res_cat_df = pd.DataFrame(index=cat_df.index, columns=cat_df.columns)

        samples = cat_df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = cat_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:
                plates = (
                    cat_df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )

                for plate in plates:

                    for idx in cat_df.index:

                        tgt_loc = (sample, tp, plate)

                        if cat_df.loc[idx, tgt_loc] == "I":
                            res_cat_df.loc[idx, tgt_loc] = 0
                        elif cat_df.loc[idx, tgt_loc] == "B":
                            res_cat_df.loc[idx, tgt_loc] = 1
                        elif cat_df.loc[idx, tgt_loc] == "P":
                            res_cat_df.loc[idx, tgt_loc] = 2
                        else:
                            res_cat_df.loc[idx, tgt_loc] == "E"

        return res_cat_df

    # * Generate dataframe containing the average category
    def generate_category_dfs(self, categories_dataframe: pd.DataFrame):

        cat_df = categories_dataframe

        res_ave_cat_df = pd.DataFrame(index=cat_df.index)

        # Iteration
        samples = cat_df.columns.get_level_values(level=0).unique()
        for sample in samples:

            tps = cat_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:

                ## Determine the mean across every sample/timepoint combination
                tgt_means = cat_df.loc[:, (sample, tp)].mean(axis=1)

                ## Category determined by the mayority, this math only works if there are 3 plates
                res_cats = []
                for mean in tgt_means:
                    if mean < 0.5:
                        res_cats.append(0)  # 0 is equal to "Indifferent"
                    elif mean < 1.5:
                        res_cats.append(1)  # 1 is equal to "Boundary"
                    else:
                        res_cats.append(2)  # 2 is equal to "Positive"

                ## Add resulting category to Results Average Categories df
                res_ave_cat_df[(sample, tp)] = res_cats

                multi_col = pd.MultiIndex.from_tuples(
                    res_ave_cat_df.columns, names=["Sample", "TP"]
                )

                res_ave_cat_df.columns = multi_col

                res_ave_cat_df = res_ave_cat_df.sort_index(axis=1)

        return res_ave_cat_df

    # * Generate dataframes containing the average category and range stats
    def generate_category_and_range_dfs(
        self,
        data_dataframe: pd.DataFrame,
        categories_dataframe: pd.DataFrame,
    ):

        # Rename df objects with shorter named variables
        data_df = data_dataframe
        cat_df = categories_dataframe

        # Copy data dataframe for its layout
        mean_cat_df = data_df.copy()

        # Generate lists to store timepoints and range statistics
        sample_tps = []  # Store (sample, timepoint) combinations
        range_mins = []  # Store minimum value across plates
        range_maxs = []  # Story maximum value across plates
        range_difs = []  # Store maximum range

        # Iteration
        samples = cat_df.columns.get_level_values(level=0).unique()
        for sample in samples:

            tps = cat_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:

                ## Determine the mean across every sample/timepoint combination
                tgt_means = cat_df.loc[:, (sample, tp)].mean(axis=1)

                ## Category determined by the mayority, this math only works if there are 3 plates
                res_cats = []
                for mean in tgt_means:
                    if mean < 0.5:
                        res_cats.append(0)  # 0 is equal to "Indifferent"
                    elif mean < 1.5:
                        res_cats.append(1)  # 1 is equal to "Boundary"
                    else:
                        res_cats.append(2)  # 2 is equal to "Positive"

                ## Column name of the established categorization
                column_title = f"Cat. Mean {str(tp)}"

                ## Add resulting category to general mean target dataframe
                mean_cat_df[(sample, tp, column_title)] = res_cats
                mean_cat_df = mean_cat_df.sort_index(axis=1)

                ## Determine the minimum and maximum value across all plates
                xy_min = data_df.loc[:, (sample, tp)].min().min()
                xy_max = data_df.loc[:, (sample, tp)].max().max()
                xy_dif = xy_max - xy_min

                ## Round up values and change to float
                smpl_tp = (sample, tp)
                xy_min = float(round(xy_min, 3))
                xy_max = float(round(xy_max, 3))
                xy_dif = float(round(xy_dif, 3))

                ## Append range stats to their respective lists
                sample_tps.append(smpl_tp)
                range_mins.append(xy_min)
                range_maxs.append(xy_max)
                range_difs.append(xy_dif)

        multi_strains_tp = pd.MultiIndex.from_tuples(sample_tps, names=["Sample", "TP"])

        # Establish dictionary for range stats to generate dataframe
        range_dict = {
            "Range Min": range_mins,
            "Range Max": range_maxs,
            "Range Dif": range_difs,
        }

        # Range Stats Dataframe
        range_stats_df = pd.DataFrame(range_dict, index=multi_strains_tp)

        return mean_cat_df, range_stats_df

    def determine_categories_df(
        self,
        categories_dataframe: pd.DataFrame,
    ):

        trans_cat_df = self.transform_ids_to_num(
            categories_dataframe=categories_dataframe,
        )

        res_cat_df = self.generate_category_dfs(
            categories_dataframe=trans_cat_df,
        )

        return res_cat_df

    def determine_categories_and_range_stats(
        self,
        data_dataframe: pd.DataFrame,
        categories_dataframe: pd.DataFrame,
    ):

        trans_cat_df = self.transform_ids_to_num(categories_dataframe)

        res_mean_cat_df, res_rng_stats_df = self.generate_category_and_range_dfs(
            data_dataframe=data_dataframe,
            categories_dataframe=trans_cat_df,
        )

        return res_mean_cat_df, res_rng_stats_df

    def values_distribution(
        self,
        data_dataframe: pd.DataFrame,
        categories_dataframe: pd.DataFrame,
        strain_id_dataframe: pd.DataFrame,
        kde_thresh: float = 0.1,
        kde_levels: int = 5,
        diag_cut: float = 2,
        save_figures: bool = False,
    ):

        # Transform categories to numbers and generate mean_categories_dataframe and range_dataframe
        data_df, range_stats_df = self.determine_categories_and_range_stats(
            data_dataframe=data_dataframe,
            categories_dataframe=categories_dataframe,
        )

        # * Series of inputs to name files if saved
        ## Input for data type
        if save_figures is True:
            data_name = input(
                "Name of the data dataframe? sample_tp_[input]_categorizationname"
            )

            if data_name == "":
                save_figures = False

        ## Input for categorization type
        if save_figures is True:
            id_name = input(
                "Name of the categorization dataframe? sample_tp_dataname_[input]"
            )

            if id_name == "":
                save_figures = False

        # * Series of inputs for figures

        ## Input for data type
        data_type = input("What is the data's type? (BWA, 535nm, 590nm, etc...)")

        ## Input for whether the data has been normalized by BNT
        norm_data = input("Is the data normalized? (Y, N)")
        if norm_data == "Y":
            norm_data = "normalized"
        else:
            norm_data = "raw"

        ## Input for categorization type
        cat_type = input("What is the categorization type (BIOLOG, BNT, BEAT)")

        # Iteration
        samples = data_df.columns.get_level_values(level=0).unique()
        for sample in samples:

            # Fetch strain name to substitute ID number
            strain_name = strain_id_dataframe.loc[sample, "Strain"]

            tps = data_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:

                # Fetch maximum and minimum from range stats df
                xy_min = range_stats_df.loc[(sample, tp), "Range Min"]
                xy_max = range_stats_df.loc[(sample, tp), "Range Max"]

                # Determine the standard maximum and minimum values to use in graphs
                xy_min = xy_min - (0.2 * xy_max)
                xy_max = xy_max + (0.2 * xy_max)

                # Set current target section of the general mean target dataframe for plotting
                tgt_in_df = data_df.loc[:, (sample, tp)]

                # Column name of the established categorization
                title = f"Cat. Mean {str(tp)}"

                # Function to plot regression line plot
                def plot_regplot(x, y, **kwargs):

                    if kwargs["label"] == 0:
                        sns.regplot(
                            data=kwargs["data"],
                            x=x.name,
                            y=y.name,
                            scatter=False,
                            color=kwargs["color"],
                            truncate=False,
                        )

                # Results are displayed in a Pair Grid system of seaborn
                g = sns.PairGrid(
                    data=tgt_in_df,
                    hue=title,
                    corner=False,
                    palette="colorblind",
                    height=4,
                )

                # Name PairGrids by strain name and timepoint
                g.figure.suptitle(
                    f"{strain_name} at {tp} hrs\nData: {data_type} ({norm_data})   Categorization: {cat_type} ",
                    y=1.02,
                    fontweight="bold",
                )

                # Diagonal Histograms
                g.map_diag(
                    sns.histplot,
                    multiple="stack",
                    element="step",
                    kde=False,
                    bins=25,
                )

                # Lower Diagonal Scatter and KDE Plots
                g.map_lower(
                    sns.scatterplot,
                )

                g.map_lower(plot_regplot, color="crimson", data=tgt_in_df)

                # Upper KDE Plots
                g.map_upper(
                    sns.kdeplot,
                    levels=kde_levels,
                    thresh=kde_thresh,
                    fill=True,
                    warn_singular=False,
                )

                # Establish minimum and maximum axis values
                g.set(
                    ylim=(xy_min, xy_max),
                    xlim=(xy_min, xy_max),
                )

                # Edit Legend Properties
                g.add_legend()

                ## Legend Title
                leg_title = "Categories"
                g._legend.set_title(leg_title)

                ## Legend Labels
                leg_labels = ["Indifferent", "Boundary", "Positive"]
                for text, label in zip(g._legend.texts, leg_labels):
                    text.set_text(label)

                # Calculate Slope, Intercept, R-Value, P-Value and Standard Error for the regression line
                axs = g.figure.get_axes()
                ## Iteration
                for ax in axs:
                    lines = ax.get_lines()
                    for line in lines:
                        if len(line.get_xydata()) > 0:

                            xd = line.get_xdata()
                            yd = line.get_ydata()

                            slope, intercept, r, p, sterr = stats.linregress(
                                x=xd,
                                y=yd,
                            )

                            # Add regression equation to plot
                            if intercept < 0:
                                formula = f"y = {str(round(slope, 2))}x - {str(round(abs(intercept),2))}"
                            else:
                                formula = f"y = {str(round(slope, 2))}x + {str(round(abs(intercept),2))}"

                            # Set equation coordinates in plot
                            ax.text(
                                x=0.05 * xy_max,
                                y=0.9 * xy_max,
                                s=formula,
                            )

                # Save figure
                if save_figures is True:

                    fig_name = (
                        f"{str(sample)}_{str(tp)}_{str(data_name)}_{str(id_name)}.png"
                    )

                    today_data_fdr = self.generate_daily_folder()
                    fig_path = today_data_fdr / fig_name
                    fig_path.resolve()

                    g.savefig(fig_path)

                plt.show()

    def determine_essential_metabolites(
        self,
        mets_dataframe: pd.DataFrame,
    ):

        mets_df = mets_dataframe

        e_mets = [[met] for met in mets_df.loc[mets_df["Essential"] == 1, "Compound"]]

        return e_mets

    # Generate a dataframe which defines range stats and metabolites of interest (positive, boundary) based on average categorization per strain/tp combination
    def generate_summary_report(
        self,
        data_dataframe: pd.DataFrame,
        categories_dataframe: pd.DataFrame,
        strain_id_dataframe: pd.DataFrame,
        mets_dataframe: pd.DataFrame,
        save_figures=False,
    ):

        # * Series of inputs to name files if saved
        ## Input for data type
        if save_figures is True:
            data_name = input(
                "Name of the data dataframe? sample_tp_[input]_categorizationname"
            )

            if data_name == "":
                save_figures = False

        ## Input for categorization type
        if save_figures is True:
            id_name = input(
                "Name of the categorization dataframe? sample_tp_dataname_[input]"
            )

            if id_name == "":
                save_figures = False

        # * Series of inputs for figures

        ## Input for data type
        data_type = input("What is the data's type? (BWA, 535nm, 590nm, etc...)")

        ## Input for categorization type
        cat_type = input("What is the categorization type (BIOLOG, BNT, BEAT)")

        # Transform categories to numbers and generate mean_categories_dataframe and range_dataframe
        mean_cat_df, range_df = self.determine_categories_and_range_stats(
            data_dataframe=data_dataframe,
            categories_dataframe=categories_dataframe,
        )

        e_mets_list = self.determine_essential_metabolites(
            mets_dataframe=mets_dataframe,
        )

        # Compare lists function
        def comp(rlist, clist):
            colors = []
            for cval in clist:
                if cval in rlist:
                    colors.append(["tab:green"])
                else:
                    colors.append(["w"])
            return colors

        # Iteration
        samples = mean_cat_df.columns.get_level_values(level=0).unique()
        for sample in samples:

            # Fetch strain name to substitute ID number
            strain_name = strain_id_dataframe.loc[sample, "Strain"]

            tps = mean_cat_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:

                title = f"Cat. Mean {str(tp)}"

                positive = [
                    [met]
                    for met in mean_cat_df.index
                    if mean_cat_df.loc[met, (sample, tp, title)] == 2
                ]
                boundary = [
                    [met]
                    for met in mean_cat_df.index
                    if mean_cat_df.loc[met, (sample, tp, title)] == 1
                ]
                indifferent = [
                    [met]
                    for met in mean_cat_df.index
                    if mean_cat_df.loc[met, (sample, tp, title)] == 0
                ]

                pie_dict = {
                    "Sample - TP": (sample, tp),
                    "Strain": strain_name,
                    "Positive": positive,
                    "Boundary": boundary,
                    "Indifferent": indifferent,
                }

                ## Creating Datasets
                categories = ["Positive", "Boundary", "Indifferent"]
                data = [len(positive), len(boundary), len(indifferent)]

                ## Creating explode data
                explode = (0.1, 0.1, 0.1)

                ## Creating color parameters
                colors = ("tab:blue", "tab:orange", "tab:gray")

                ## Wedge properties
                wp = {
                    "linewidth": 1,
                    "edgecolor": "black",
                }

                def autocpt(pct, allvalues):
                    absolute = int(round(pct * np.sum(allvalues) / 100))
                    return "{:.1f}%\n({:d})".format(pct, absolute)

                fig = plt.figure(figsize=(16, 12))
                fig.set_layout_engine("tight")
                ax_1 = fig.add_subplot(1, 4, 1)

                wedges, texts, autotexts = ax_1.pie(
                    data,
                    autopct=lambda pct: autocpt(pct, data),
                    explode=explode,
                    shadow=False,
                    colors=colors,
                    startangle=90,
                    wedgeprops=wp,
                    textprops=dict(color="black"),
                )

                ## Adding legend
                ax_1.legend(
                    wedges,
                    categories,
                    title="Categories",
                    loc="center right",
                    bbox_to_anchor=(0.1, -0.6, 0.5, 1),
                )

                plt.setp(
                    autotexts,
                    size=12,
                    weight="bold",
                )

                ## Assign colors for positive and boundary cells if metabolites present in essential metabolites list
                p_colors = comp(rlist=e_mets_list, clist=positive)
                b_colors = comp(rlist=e_mets_list, clist=boundary)

                # Table essential metabolites
                if len(e_mets_list) != 0:
                    ax_2 = fig.add_subplot(1, 4, 2)
                    ax_2.table(
                        cellText=e_mets_list,
                        colLabels=["Essential"],
                        colColours=["tab:green"],
                        loc="upper center",
                        cellLoc="center",
                        edges="closed",
                    ).scale(0.9, 1.3)
                    ax_2.axis("off")

                # Table positive metabolites
                if len(positive) != 0:
                    ax_3 = fig.add_subplot(1, 4, 3)
                    ax_3.table(
                        cellText=positive,
                        colLabels=["Positive"],
                        colColours=["tab:blue"],
                        cellColours=p_colors,
                        loc="upper center",
                        cellLoc="center",
                        edges="closed",
                    ).scale(0.9, 1.3)
                    ax_3.axis("off")

                # Table boundary metabolites
                if len(boundary) != 0:
                    ax_4 = fig.add_subplot(1, 4, 4)
                    ax_4.table(
                        cellText=boundary,
                        colLabels=["Boundary"],
                        colColours=["tab:orange"],
                        cellColours=b_colors,
                        loc="upper center",
                        cellLoc="center",
                        edges="closed",
                    ).scale(0.9, 1.3)
                    ax_4.axis("off")

                fig.suptitle(
                    f"{strain_name} at {tp} hrs\nData: {data_type}   Categorization: {cat_type}",
                    y=0.97,
                    fontweight="bold",
                )

                # Save figure
                if save_figures is True:

                    fig_name = f"sr_{str(sample)}_{str(tp)}_{str(data_name)}_{str(id_name)}.png"

                    today_data_fdr = self.generate_daily_folder()
                    fig_path = today_data_fdr / fig_name
                    fig_path.resolve()

                    fig.savefig(fig_path)

                plt.show()

    # Generate a new folder named by Year_Month_Day if not available, get its path
    def generate_daily_folder(self) -> Path:

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

    # * GENERATE SUPER REPORT
    def generate_super_report(
        self,
        cat_biolog_df: pd.DataFrame,
        cat_bwa_norm_df: pd.DataFrame,
        cat_590nm_norm_df: pd.DataFrame,
        strain_id_df: pd.DataFrame,
        met_assay_df: pd.DataFrame,
        width: float = 6,
        height: float = 20,
        save_figures: bool = False,
    ):

        # Generate Average Categorization for BIOLOG Categorizations
        ave_cat_biolog_df = self.determine_categories_df(
            categories_dataframe=cat_biolog_df,
        )

        # Generate Average Categorization for BWA Normalized Categorizations
        ave_cat_bwa_df = self.determine_categories_df(
            categories_dataframe=cat_bwa_norm_df,
        )

        # Generate Average Categorization for 590nm Normalized Categorizations
        ave_cat_590nm_df = self.determine_categories_df(
            categories_dataframe=cat_590nm_norm_df,
        )

        def assign_color(tgt_list):
            colors = []
            for tgt in tgt_list:
                if tgt == 2:
                    colors.append("tab:blue")
                elif tgt == 1:
                    colors.append("tab:orange")
                else:
                    colors.append("w")

            return colors

        def assign_met_assay_color(tgt_list):
            colors = []
            for tgt in tgt_list:
                if tgt >= 0.2:
                    colors.append("tab:blue")
                elif tgt < 0.2 and tgt >= 0.1:
                    colors.append("tab:orange")
                else:
                    colors.append("w")

            return colors

        def combine_colors(
            c_biolog,
            c_bwa,
            c_590nm,
            c_fba,
        ):

            n = len(c_biolog)
            c_comb = []

            for i in range(n):
                i_c = []

                i_c.append(c_biolog[i])
                i_c.append(c_bwa[i])
                i_c.append(c_590nm[i])
                i_c.append(c_fba[i])

                c_comb.append(i_c)

            return c_comb

        # Iteration

        samples = ave_cat_biolog_df.columns.get_level_values(level=0).unique()
        for sample in samples:

            # summary_df = pd.DataFrame(index=ave_cat_biolog_df.index)

            # Fetch strain name to substitute ID number
            strain_name = strain_id_df.loc[sample, "Strain"]
            strain_dsmz = strain_id_df.loc[sample, "DSMZ-number"]

            tgt_mAssay = met_assay_df.loc[:, strain_dsmz]

            tps = (
                ave_cat_biolog_df.loc[:, sample]
                .columns.get_level_values(level=0)
                .unique()
            )
            for tp in tps:

                tgt_biolog = ave_cat_biolog_df.loc[:, (sample, tp)]
                tgt_bwa = ave_cat_bwa_df.loc[:, (sample, tp)]
                tgt_590nm = ave_cat_590nm_df.loc[:, (sample, tp)]

                summary_df = pd.DataFrame(index=ave_cat_biolog_df.index)
                summary_df[("BIOLOG")] = tgt_biolog
                summary_df[("BWA")] = tgt_bwa
                summary_df[("590nm")] = tgt_590nm
                summary_df[("FBA")] = tgt_mAssay

                summary_df = summary_df.sort_values(
                    summary_df.columns.tolist(),
                    ascending=False,
                )

                colors_biolog = assign_color(summary_df["BIOLOG"])
                colors_bwa = assign_color(summary_df["BWA"])
                colors_590nm = assign_color(summary_df["590nm"])
                colors_fba = assign_met_assay_color(summary_df["FBA"])

                colors = combine_colors(
                    c_biolog=colors_biolog,
                    c_bwa=colors_bwa,
                    c_590nm=colors_590nm,
                    c_fba=colors_fba,
                )

                fig, ax = plt.subplots(figsize=(width, height))

                table1 = ax.table(
                    cellText=summary_df.values,
                    cellColours=colors,
                    rowLabels=summary_df.index,
                    rowLoc="right",
                    colLabels=summary_df.columns,
                    cellLoc="center",
                    loc="best",
                )

                table1.auto_set_font_size(False)
                table1.set_fontsize(12)

                fig.suptitle(f"{strain_name} {tp}")
                ax.axis("off")
                fig.tight_layout()

                # Save figure
                if save_figures is True:

                    fig_name = f"{str(sample)}_{str(tp)}_sumtable.png"

                    today_data_fdr = self.generate_daily_folder()
                    fig_path = today_data_fdr / fig_name
                    fig_path.resolve()

                    fig.savefig(fig_path)

                plt.show()
