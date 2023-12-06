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

    def yeojohn_normalization(
        self,
    ) -> (
        pd.DataFrame
    ):  # function to normalize the difference dataframe with yeojohnson transformation approach
        self.calculation_dif()

        self.norm_dict = dict()

        yeojohnTr = PowerTransformer(
            standardize=True
        )  # yeo-johnson transformation is the default method attribute

        for column in self.dif_df.columns:
            norm_data = yeojohnTr.fit_transform(
                self.dif_df[column].values.reshape(-1, 1)
            )

            norm_data = norm_data.reshape(1, len(norm_data))[0]

            self.norm_dict[column] = norm_data

        self.normalized_df = pd.DataFrame.from_dict(self.norm_dict)

        self.normalized_df.index = self.dif_df.index

        return self.normalized_df

    def generate_skewness_data(self) -> pd.DataFrame:
        self.yeojohn_normalization()
        old_skew = self.dif_df.skew()
        new_skew = self.normalized_df.skew()

        skew_dict = {"Before Trans.": old_skew, "After Trans.": new_skew}

        self.skew_df = pd.DataFrame.from_dict(skew_dict)
        self.skew_df.Name = "Skewness Information"

        return self.skew_df

    def style_boundary(self, value):
        std = 1.005249

        boundary = [1.5 * std < v <= 2 * std for v in value]

        return [f"color: tomato" if i else None for i in boundary]

    def style_positive(self, value):
        std = 1.005249

        positive = value > 2 * std

        return [f"color: royalblue" if i else None for i in positive]

    def display_data_categorization(self, dataframe):
        result = dataframe.style.apply(self.style_positive, axis=1).apply(
            self.style_boundary, axis=1
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
