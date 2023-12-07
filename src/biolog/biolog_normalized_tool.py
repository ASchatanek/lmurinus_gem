import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from sklearn.preprocessing import PowerTransformer

sns.set_theme()
sns.set_palette(palette="rainbow")


class BNT:
    def __init__(self, dataframe: pd.DataFrame) -> None:
        self.df = dataframe

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

    def yeojohn_normalization(
        self,
    ) -> (
        pd.DataFrame
    ):  # function to normalize the difference dataframe with yeojohnson transformation approach
        self.calculation_dif()
        self.calculation_ave()
        self.add_ave_dataframe()

        self.norm_dict = dict()

        yeojohnTr = PowerTransformer(
            standardize=True
        )  # yeo-johnson transformation is the default method attribute

        for column in self.complete_df.columns:
            norm_data = yeojohnTr.fit_transform(
                self.complete_df[column].values.reshape(-1, 1)
            )

            norm_data = norm_data.reshape(1, len(norm_data))[0]

            self.norm_dict[column] = norm_data

        self.normalized_df = pd.DataFrame.from_dict(self.norm_dict)

        self.normalized_df.index = self.complete_df.index

        return self.normalized_df

    def generate_skewness_data(self) -> pd.DataFrame:
        self.yeojohn_normalization()
        old_skew = self.dif_df.skew()
        new_skew = self.normalized_df.skew()

        skew_dict = {"Before Trans.": old_skew, "After Trans.": new_skew}

        self.skew_df = pd.DataFrame.from_dict(skew_dict)
        self.skew_df.Name = "Skewness Information"

        return self.skew_df

    def style_categories(self, column, limit_i=1.5, limit_b=2):
        std = 1.005249
        result = []

        for cell in column:
            if cell <= limit_i * std:
                result.append("background-color: white;color: darkgrey")
            elif cell > limit_i and cell <= limit_b:
                result.append(
                    "background-color: khaki; font-weight: bold; color: goldenrod"
                )
            else:
                result.append(
                    "background-color: palegreen; font-weight: bold; color: forestgreen"
                )

        return result

    def display_categories(self, dataframe: pd.DataFrame, limit_i=1.5, limit_b=2):
        result = dataframe.style.set_table_attributes("style='display:inline'").apply(
            self.style_categories, limit_i=limit_i, limit_b=limit_b, axis=1
        )

        return result

    def display_yeojohn_norm(self, target_columns=None):
        self.yeojohn_normalization()

        if target_columns != None:
            target_norm_df = self.normalized_df[target_columns]
            target_dif_df = self.dif_df[target_columns]

            for column in target_norm_df:
                plt.figure(figsize=(30, 6))

                plt.subplot(1, 4, 1)
                plt.title("Distribution before Transformation", fontsize=15)
                sns.histplot(target_dif_df[column], kde=True, color="red")

                plt.subplot(1, 4, 2)
                stats.probplot(target_dif_df[column], dist="norm", plot=plt)  # QQ Plot
                plt.title("QQ Plot before Transformation", fontsize=15)

                plt.subplot(1, 4, 3)
                plt.title("Distribution after Transformation", fontsize=15)
                sns.histplot(target_norm_df[column], bins=20, kde=True, legend=False)

                plt.subplot(1, 4, 4)
                stats.probplot(target_norm_df[column], dist="norm", plot=plt)  # QQ Plot
                plt.title("QQ Plot after Transformation", fontsize=15)

                plt.show()

        else:
            for column in self.normalized_df.columns:
                plt.figure(figsize=(30, 6))

                plt.subplot(1, 4, 1)
                plt.title("Distribution before Transformation", fontsize=15)
                sns.histplot(self.dif_df[column], kde=True, color="red")

                plt.subplot(1, 4, 2)
                stats.probplot(self.dif_df[column], dist="norm", plot=plt)  # QQ Plot
                plt.title("QQ Plot before Transformation", fontsize=15)

                plt.subplot(1, 4, 3)
                plt.title("Distribution after Transformation", fontsize=15)
                sns.histplot(
                    self.normalized_df[column], bins=20, kde=True, legend=False
                )

                plt.subplot(1, 4, 4)
                stats.probplot(
                    self.normalized_df[column], dist="norm", plot=plt
                )  # QQ Plot
                plt.title("QQ Plot after Transformation", fontsize=15)

                plt.show()
