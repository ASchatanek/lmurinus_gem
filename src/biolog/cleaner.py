from pathlib import Path
import pandas as pd


class Cleaner:
    def rinse(self):
        self.core = []
        self.index = []
        self.col_list = list(self.df.columns.values)

        for col in self.col_list:
            if len(col) > 3:
                self.index.append(col)
            elif len(col) == 3:
                self.core.append(col)

        if len(self.core) == len(self.compounds):
            self.new_cols = self.index + self.compounds

    def organize(self):
        self.df.columns = self.new_cols
        self.df.set_index(["Sample #", "Plate #", "Hours"], inplace=True)
        self.df.index.names = ["Sample #", "Plate #", "TP"]
        self.df = (
            self.df.transpose()
            .reorder_levels(["Sample #", "TP", "Plate #"], axis=1)
            .sort_index(axis=1)
        )

    def cleaNorganize(self, dataframe, compounds) -> pd.DataFrame:
        self.df = dataframe
        self.compounds = list(compounds)
        self.rinse()
        self.organize()
        return self.df
