import pandas as pd
import numpy as np
import matplotlib.pylab as plt


class BAT:
    def __init__(self, df):
        self.df = df

    def data_analysis(self):
        self.calculation_dif()
        self.calculation_ave()
        self.calculation_bounds()
        self.calculation_stats()
        self.add_ave_dataframe()
        self.identify_result_type()

    def calculation_dif(self) -> pd.DataFrame:
        # Define variables for the process
        self.dif_df = pd.DataFrame(index=self.df.index, columns=self.df.columns)

        # Use the information of the dataset hierarchy to calculate the difference of each value to the "Water" value of each column
        samples = self.df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = self.df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:
                plates = (
                    self.df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )
                for plate in plates:
                    dif = (
                        self.df.loc[:, (sample, tp, plate)]
                        - self.df.loc["Water", (sample, tp, plate)]
                    )
                    self.dif_df[(sample, tp, plate)] = dif

        self.dif_df = self.dif_df.sort_index(axis=1)

        return self.dif_df

    def calculation_ave(self) -> pd.DataFrame:
        # Define variables for the process
        self.tp_tuples = []

        # Use the information of the hierarchy system to define the tuple list and dictionary
        samples = self.df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = self.df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:
                self.tp_tuples.append((sample, tp))

        # Create MultiIndex columns
        self.multi_columns = pd.MultiIndex.from_tuples(
            self.tp_tuples, names=["Sample #", "TP"]
        )

        # Create empty dataframe with correct index and columns
        self.ave_df = pd.DataFrame(index=self.df.index, columns=self.multi_columns)

        # Calculate the average of plates for each row in the dataset
        samples = self.ave_df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = self.ave_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:
                tp_average = self.dif_df.loc[:, (sample, tp)].mean(axis=1)
                self.ave_df[sample, tp] = tp_average

        # Sort and round columns of resulting dataframe
        self.ave_df = self.ave_df.sort_index(axis=1)
        # self.ave_df = self.ave_df.round(2)

        # Returns dataframe with the average of the plates for each timepoint of each sample
        return self.ave_df

    def calculation_bounds(self) -> pd.DataFrame:
        # Get data description, transpose columns to index and sort index
        self.des_df = self.ave_df.describe().transpose().sort_index()

        # Calculate lower and upper bounds for each dataset
        self.percentages = ["25%", "50%", "75%"]
        self.df_info = self.des_df[self.percentages].copy()

        ## Calculate IQR and add IQR column
        IQR = (self.df_info.loc[:, ("75%")] - self.df_info.loc[:, ("25%")]).abs()
        self.df_info.loc[:, ("IQR")] = IQR

        ## Calculate upper bound (UB) and lower bound (LB), add columns
        UB = self.df_info.loc[:, ("75%")] + (self.df_info.loc[:, ("IQR")] * 1.5)
        self.df_info.loc[:, ("UB")] = UB

        LB = self.df_info.loc[:, ("25%")] - (self.df_info.loc[:, ("IQR")] * 1.5)
        self.df_info.loc[:, ("LB")] = LB

        return self.df_info

    def calculation_stats(self) -> pd.DataFrame:
        # * The idea is that a dataframe containing the B&W stats, IQR, boundaries and the std. devs. of each plate and their average.

        # Define the variables using defined functions for average and boundary calculation
        self.bds_df = self.calculation_bounds()

        # Modifying boundaries dataframe's format
        self.bds_df = (
            self.bds_df.transpose()
        )  # Transpose the structure of the boundary dataframe and add the name "Avg. Stats."

        # Define tuples and categories to define the new multiindex structure
        self.idx_tups = []
        self.idx_cat = ["B&W Stats", "Average", "Std. Dev."]

        # Create a tuple that contains the numbering of each plate
        #! This assumes there are only 3 plates, this will need to change if the amount of plates varies
        for i in range(1, 4):
            self.idx_tups.append((self.idx_cat[1], f"Plate {i}"))
        self.idx_tups.append((self.idx_cat[1], "Plate Avg."))

        for i in range(1, 4):
            self.idx_tups.append((self.idx_cat[2], f"Plate {i}"))
        self.idx_tups.append((self.idx_cat[2], "Plate Avg."))

        # Define multiindex with generated list of tuples
        self.idx_mult = pd.MultiIndex.from_tuples(self.idx_tups)

        # Create a concatenated dataframe from the newly designed stats_df dataframe and existing bds_df (containing the calculated B&W statistics)
        ## Concatenation
        self.stats_df = pd.DataFrame(index=self.idx_mult, columns=self.bds_df.columns)
        self.stats_df = pd.concat([self.bds_df, self.stats_df])

        ## Rename the boundary index names to a tuple containing the levels for a multiindex
        for value in self.stats_df.index.values:
            if type(value) == str:
                self.stats_df = self.stats_df.rename(
                    index={value: (self.idx_cat[0], value)}
                )

        ## Take the generated tuples to create the new multiindex
        self.stats_df.index = pd.MultiIndex.from_tuples(self.stats_df.index.values)

        # * The Std. Dev. of each plate is calculated. This is done after filtering the values within each plate based on the lower and upper bound defined of their average B&W whisker´s range.
        #! The database must follow a multiindex system. They should be organized as follows: "Sample #" --> "TP" --> "Plate #"
        self.df_samples = self.df.columns.get_level_values(
            level=0
        ).unique()  # Define the samples in the desired database

        ## Iterate each sample in the database
        for sample in self.df_samples:
            tps = self.dif_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:
                ub = self.stats_df.loc[("B&W Stats", "UB"), (sample, tp)]  # Upper Bound
                lb = self.stats_df.loc[("B&W Stats", "LB"), (sample, tp)]  # Lower Bound

                plates = (
                    self.dif_df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )

                # Iterate each plate, filter values based on the defined upper and lower bounds and calculate their std. dev.
                for plate in plates:
                    ## Identify target
                    target = self.dif_df.loc[:, (sample, tp, plate)]
                    ## Assign result by filtering for values within the defined bounds
                    result = target[(lb < target) & (target < ub)]
                    ## Calculate filtered values' mean and std. dev.
                    self.stats_df.loc[("Average", plate), (sample, tp)] = result.mean()
                    self.stats_df.loc[("Std. Dev.", plate), (sample, tp)] = result.std()

                # Once each individuals plates' std. dev. is calculated, the std. dev. of values average is calculated (after filtering)
                ## Identify target
                target = self.ave_df.loc[:, (sample, tp)]
                ## Assign result by filtering for values within the defined bounds
                result = target[(lb < target) & (target < ub)]
                ## Calculate filtered values' mean and std. dev.
                self.stats_df.loc[
                    ("Average", "Plate Avg."), (sample, tp)
                ] = result.mean()
                self.stats_df.loc[
                    ("Std. Dev.", "Plate Avg."), (sample, tp)
                ] = result.std()

        return self.stats_df

    def add_ave_dataframe(self) -> pd.DataFrame:
        self.complete_df = self.dif_df.copy()

        samples = self.complete_df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = (
                self.complete_df.loc[:, sample]
                .columns.get_level_values(level=0)
                .unique()
            )
            for tp in tps:
                self.complete_df.loc[:, (sample, tp, "Plate Avg.")] = self.ave_df.loc[
                    :, (sample, tp)
                ]

        self.complete_df = self.complete_df.sort_index(axis=1)

        return self.complete_df

    def identify_result_type(self) -> pd.DataFrame:
        tags = ["I", "B", "P"]

        self.complete_df = self.add_ave_dataframe()

        self.id_df = pd.DataFrame(
            index=self.complete_df.index, columns=self.complete_df.columns
        )

        samples = self.complete_df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = self.id_df.loc[:, sample].columns.get_level_values(level=0).unique()
            for tp in tps:
                plates = (
                    self.id_df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )
                for plate in plates:
                    target = self.complete_df.loc[:, (sample, tp, plate)]
                    target_mean = self.stats_df.loc[("Average", plate), (sample, tp)]
                    target_stdev = self.stats_df.loc[("Std. Dev.", plate), (sample, tp)]

                    limit_indi = target_mean + (target_stdev * 2)
                    limit_bdry = target_mean + (target_stdev * 4)

                    self.stats_df.loc[
                        ("Limit Indiferent", plate), (sample, tp)
                    ] = limit_indi

                    self.stats_df.loc[
                        ("Limit Boundary", plate), (sample, tp)
                    ] = limit_bdry

                    conditions = [
                        (target <= limit_indi),
                        (target > limit_indi) & (target <= limit_bdry),
                        (target > limit_bdry),
                    ]

                    self.id_df[(sample, tp, plate)] = np.select(conditions, tags)

        stats_df_idx_order = [
            "B&W Stats",
            "Average",
            "Std. Dev.",
            "Limit Indiferent",
            "Limit Boundary",
        ]
        self.stats_df = self.stats_df.sort_index().reindex(stats_df_idx_order, level=0)

        return self.id_df

    def add_tags(self, value):
        if value == "I":
            return "color : darkgrey"
        elif value == "B":
            return "color : goldenrod"  # darkred
        elif value == "P":
            return "color : forestgreen"  # darkblue

    def add_highlight(self, value):
        if value == "I":
            return "background-color : white"
        elif value == "B":
            return "background-color : khaki"  # lightcoral
        elif value == "P":
            return "background-color : palegreen"  # lightblue

    def bold_values(self, value):
        if value == "I":
            return "color : darkgrey"
        elif value == "B":
            return "font-weight: bold"
        elif value == "P":
            return "font-weight: bold"

    def style_categories(self, value):
        if value == "I":
            return "color: darkgrey"
        elif value == "B":
            return "background-color: khaki; font-weight: bold; color: goldenrod"
        elif value == "P":
            return "background-color: palegreen; font-weight: bold; color: forestgreen"

    def display_categories(self, dataframe: pd.DataFrame, tag_df=None, columns=None):
        if tag_df is None:
            tag_df = self.id_df

        if columns is None:
            dataframe = dataframe
            return (
                dataframe.style.set_properties(**{"background-color": "white"})
                .apply(lambda x: tag_df.map(self.style_categories), axis=None)
                .format(precision=2)
            )
        else:
            dataframe = dataframe.loc[:, columns]
            return (
                dataframe.style.set_properties(**{"background-color": "white"})
                .apply(
                    lambda x: tag_df.loc[:, columns].map(self.style_categories),
                    axis=None,
                )
                .format(precision=2)
            )

    def apply_style(self, df: pd.DataFrame, tag_df=None, columns=None):
        if tag_df is None:
            tag_df = self.id_df

        if columns is None:
            dataframe = df
            return (
                dataframe.style.set_properties(**{"background-color": "white"})
                .apply(lambda x: tag_df.applymap(self.add_tags), axis=None)
                .apply(lambda x: tag_df.applymap(self.add_highlight), axis=None)
                .apply(lambda v: tag_df.applymap(self.bold_values), axis=None)
                .format(precision=2)
            )
        else:
            dataframe = df.loc[:, columns]
            return (
                dataframe.style.set_properties(**{"background-color": "white"})
                .apply(
                    lambda x: tag_df.loc[:, columns].applymap(self.add_tags), axis=None
                )
                .apply(
                    lambda x: tag_df.loc[:, columns].applymap(self.add_highlight),
                    axis=None,
                )
                .apply(
                    lambda v: tag_df.loc[:, columns].applymap(self.bold_values),
                    axis=None,
                )
                .format(precision=2)
            )

    def set_for_display(self, df):
        target = df
        if type(target) == pd.DataFrame:
            result = target.style.set_table_attributes(
                "style='display:inline; margin-right:100px;'"
            )
            return result
        else:
            target.set_table_attributes("style='display:inline; margin-right:100px;'")
            return target

    def comparison_analysis(self, comparison_df) -> pd.DataFrame:
        self.mismatch_df = pd.DataFrame(
            columns=[
                "Target Index",
                "Compound",
                "Result -> Comparison",
                "Raw Value",
                "Plate Blank",
                "Difference",
                "Limit I",
                "Limit B",
                "Per. I",
                "Per. B",
            ]
        )

        self.comparison_df = comparison_df
        samples = self.comparison_df.columns.get_level_values(level=0).unique()
        for sample in samples:
            tps = (
                self.comparison_df.loc[:, sample]
                .columns.get_level_values(level=0)
                .unique()
            )
            for tp in tps:
                plates = (
                    self.comparison_df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )
                for plate in plates:
                    self.target_coordinates = (sample, tp, plate)

                    self.comparison_target = self.comparison_df.loc[
                        :, self.target_coordinates
                    ].str.strip()
                    self.id_df_target = self.id_df.loc[
                        :, self.target_coordinates
                    ].str.strip()

                    self.comparison = self.comparison_target != self.id_df_target

                    self.comparison_result = self.comparison_df.loc[
                        self.comparison, self.target_coordinates
                    ]

                    # Add results along others to the Mismatch Table
                    for compound in self.comparison_result.index:
                        target_data = {
                            "Target Index": self.target_coordinates,
                            "Compound": compound,
                            "Result -> Comparison": self.id_df.loc[
                                compound, self.target_coordinates
                            ]
                            + " -> "
                            + self.comparison_df.loc[compound, self.target_coordinates],
                            "Raw Value": self.df.loc[compound, self.target_coordinates],
                            "Plate Blank": self.df.loc[
                                "Water", self.target_coordinates
                            ],
                            "Difference": self.dif_df.loc[
                                compound, self.target_coordinates
                            ],
                            "Limit I": self.stats_df.loc[
                                ("Limit Indiferent", plate), (sample, tp)
                            ],
                            "Limit B": self.stats_df.loc[
                                ("Limit Boundary", plate), (sample, tp)
                            ],
                        }

                        self.mismatch_df = self.mismatch_df.append(
                            target_data, ignore_index=True
                        )

        self.mismatch_df = self.mismatch_df.set_index("Target Index")

        self.mismatch_df.index = pd.MultiIndex.from_tuples(
            self.mismatch_df.index, names=["Sample #", "Timepoint", "Plate #"]
        )

        self.limits_dis = abs(self.mismatch_df["Limit B"] - self.mismatch_df["Limit I"])

        self.mismatch_df["Per. I"] = (
            abs(self.mismatch_df["Limit I"] - self.mismatch_df["Difference"])
            / self.limits_dis
            * 100
        )
        self.mismatch_df["Per. B"] = (
            abs(self.mismatch_df["Limit B"] - self.mismatch_df["Difference"])
            / self.limits_dis
            * 100
        )

        self.mismatch_df = self.mismatch_df.sort_index(axis=0)

        ## df.loc[df[‘column’] condition, ‘new column name’] = ‘value if condition is met’
        self.mismatch_df.loc[self.mismatch_df["Per. I"] > 100.0, "Per. I"] = "-"
        self.mismatch_df.loc[self.mismatch_df["Per. B"] > 100.0, "Per. B"] = "-"

        return self.mismatch_df

    def comparison_report(self):
        self.result_target = self.mismatch_df.loc[
            :, ("Result -> Comparison", "Compound")
        ]

        self.pivot_rescomp = self.result_target.pivot_table(
            index=self.result_target.index.names,
            columns=["Result -> Comparison"],
            aggfunc="count",
            fill_value=0,
        ).droplevel(axis=1, level=0)

        self.pivot_samples = self.pivot_rescomp.index.get_level_values(level=0).unique()

        for sample in self.pivot_samples:
            tps = (
                self.pivot_rescomp.loc[sample, :]
                .index.get_level_values(level=0)
                .unique()
            )

            self.ax = self.pivot_rescomp.loc[sample, :].plot(
                kind="bar",
                title=sample,
                xlabel="Timepoints and Plates #",
                ylabel="Frequency",
                rot=30,
                figsize=(10, 10),
            )
            self.ax.tick_params(axis="x", labelsize=12)
            self.ax.tick_params(axis="y", labelsize=15)
            self.ax.title.set_size(20)

            plt.show()

            for tp in tps:
                self.series_per_I = self.mismatch_df.loc[(sample, tp), "Per. I"].copy()
                self.series_per_I = self.series_per_I[self.series_per_I != "-"].plot(
                    kind="hist",
                    bins=20,
                    xlim=(0, 100),
                    title=f"{sample, tp}",
                    xlabel="Percentages",
                )

                plt.show()

                self.pivot_sample_tp = self.result_target.loc[(sample, tp), :]

                print(sample, tp)

                self.rc_combis = self.pivot_sample_tp["Result -> Comparison"].unique()

                for comb in self.rc_combis:
                    self.pivot_s_t_comb = self.pivot_sample_tp.loc[
                        self.pivot_sample_tp["Result -> Comparison"] == comb, "Compound"
                    ]

                    print(comb)

                    self.s_t_compounds = np.unique(self.pivot_s_t_comb.values)

                    print(f"Total Mismatched Compounds: {len(self.s_t_compounds)}")

                    for count, value in enumerate(self.s_t_compounds, start=1):
                        comb_index = list(
                            self.pivot_s_t_comb[self.pivot_s_t_comb == value].index
                        )

                        if len(comb_index) == 3:
                            print(f"{count}) {value} --- All 3 Plates")
                        else:
                            print(f"{count}) {value} --- {comb_index}")

                    print("")
