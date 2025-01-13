import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy import stats


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
        id_dataframe: pd.DataFrame,
    ) -> pd.DataFrame:

        id_df = id_dataframe

        ref_id_df = pd.DataFrame(index=id_df.index, columns=id_df.columns)

        samples = id_df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = id_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:
                plates = (
                    id_df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )

                for plate in plates:

                    for idx in id_df.index:

                        tgt_loc = (sample, tp, plate)

                        if id_df.loc[idx, tgt_loc] == "I":
                            ref_id_df.loc[idx, tgt_loc] = 0
                        elif id_df.loc[idx, tgt_loc] == "B":
                            ref_id_df.loc[idx, tgt_loc] = 1
                        elif id_df.loc[idx, tgt_loc] == "P":
                            ref_id_df.loc[idx, tgt_loc] = 2
                        else:
                            ref_id_df.loc[idx, tgt_loc] == "E"

        return ref_id_df

    def compare_values_distribution(
        self,
        tgt_values_dataframe: pd.DataFrame,
        categorization_dataframe: pd.DataFrame,
        strain_id_df: pd.DataFrame,
        kde_thresh: float = 0.05,
        kde_levels: int = 10,
        save_figs: bool = False,
    ):

        tgt_df = tgt_values_dataframe
        tgt_cat_df = categorization_dataframe

        mean_tgt_df = tgt_df.copy()

        if save_figs is True:
            data_name = input(
                "Name of the data dataframe? sample_tp_[input]_categorizationname"
            )

            if data_name == "":
                save_figs = False

        if save_figs is True:
            id_name = input(
                "Name of the categorization dataframe? sample_tp_dataname_[input]"
            )

            if id_name == "":
                save_figs = False

        # Iteration
        samples = tgt_cat_df.columns.get_level_values(level=0).unique()
        for sample in samples:

            # Fetch strain name to substitute ID number
            strain_name = strain_id_df.loc[sample, "Strain"]

            tps = tgt_cat_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:

                tgt_means = tgt_cat_df.loc[:, (sample, tp)].mean(axis=1)

                res_cats = []
                for mean in tgt_means:
                    if mean < 0.5:
                        res_cats.append(0)  # Here 0 is equal to "Indifferent"
                    elif mean < 1.5:
                        res_cats.append(1)  # Here 1 is equal to "Boundary"
                    else:
                        res_cats.append(2)  # Here 2 is equal to "Positive"

                title = f"Cat. Mean {str(tp)}"

                xy_min = tgt_df.loc[:, (sample, tp)].min().min()
                xy_max = tgt_df.loc[:, (sample, tp)].max().max()

                range_dif = xy_max - xy_min

                print((sample, tp))
                print("-----------------")
                print(
                    f"min: {float(round(xy_min, 3))}  max: {float(round(xy_max, 3))}  dif: {float(round(range_dif, 3))}"
                )
                print("-----------------")

                xy_min = xy_min - (0.2 * xy_max)
                xy_max = xy_max + (0.2 * xy_max)

                mean_tgt_df[(sample, tp, title)] = res_cats
                mean_tgt_df = mean_tgt_df.sort_index(axis=1)

                tgt_in_df = mean_tgt_df.loc[:, (sample, tp)]

                def plot_regplot(x, y, **kwargs):

                    # Create Regression Line
                    if kwargs["label"] == 0:
                        sns.regplot(
                            data=kwargs["data"],
                            x=x.name,
                            y=y.name,
                            scatter=False,
                            color=kwargs["color"],
                            truncate=False,
                        )

                # Results are displayed in a pair grid system
                g = sns.PairGrid(
                    data=tgt_in_df,
                    hue=title,
                    corner=False,
                    palette="colorblind",
                    height=4,
                )

                # Name PairGrids by strain name and timepoint
                g.figure.suptitle(
                    f"{strain_name} at {tp} hrs",
                    y=1.01,
                )

                # Diagonal Histograms
                # g.map_diag(
                #     sns.histplot,
                #     multiple="stack",
                #     element="step",
                #     kde=True,
                #     bins=20,
                # )

                g.map_diag(
                    sns.kdeplot,
                    fill=True,
                    bw_adjust=2,
                    warn_singular=False,
                )

                # # Lower KDE Plots
                # g.map_lower(
                #     sns.kdeplot,
                #     levels=1,
                #     thresh=0.1,
                #     warn_singular=False,
                # )

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

                axs = g.figure.get_axes()

                # Calculate Slope, Intercept, R-Value, P-Value and Standard Error for the regression line
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

                            ax.text(
                                x=0.05 * xy_max,
                                y=0.9 * xy_max,
                                s=formula,
                            )

                # Save figure
                if save_figs is True:

                    fig_name = (
                        f"{str(sample)}_{str(tp)}_{str(data_name)}_{str(id_name)}.png"
                    )

                    today_data_fdr = self.generate_daily_folder()
                    fig_path = today_data_fdr / fig_name
                    fig_path.resolve()

                    g.savefig(fig_path)

                plt.show()

    # Function takes both changes the id dataframes to a number system and generates the pair grid distribution anaylsis figures
    def distribution_analysis(
        self,
        id_df: pd.DataFrame,
        tgt_df: pd.DataFrame,
        strain_id_df: pd.DataFrame,
        kde_thresh: float = 0.1,
        kde_levels: int = 5,
        save_figs: bool = False,
    ):

        trans_id_df = self.transform_ids_to_num(id_dataframe=id_df)

        self.compare_values_distribution(
            tgt_values_dataframe=tgt_df,
            categorization_dataframe=trans_id_df,
            strain_id_df=strain_id_df,
            kde_thresh=kde_thresh,
            kde_levels=kde_levels,
            save_figs=save_figs,
        )

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
