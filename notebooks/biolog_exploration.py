import set_cwd
import pandas_settings
import pandas as pd
import numpy as np

from src import Reader, Cleaner, BAT, BNT

reader = Reader()

biolog_values = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="BIOLOG")
abs_values = reader.find_dataframe("plate_readout_numbers.xlsx", sheet_name="Abs_535")
biolog_compounds_df = reader.find_dataframe(
    "plate_readout_numbers.xlsx", sheet_name="Compounds"
)


biolog_compounds = reader.get_compounds_list(biolog_compounds_df)

biolog_ids = reader.find_dataframe(
    "plate_readout_categorizations.xlsx", sheet_name="BIOLOG"
)

df_cleaner = Cleaner()
values_df = df_cleaner.cleaNorganize(biolog_values, compounds=biolog_compounds)
biolog_raw_id_df = df_cleaner.cleaNorganize(biolog_ids, compounds=biolog_compounds)
abs_df = df_cleaner.cleaNorganize(abs_values, compounds=biolog_compounds)

cc = BAT(values_df)
cc.data_analysis()

cc2 = BAT(abs_df)
cc2.data_analysis()

test = BNT(dataframe=values_df)
test2 = BNT(dataframe=abs_df)


df = test.yeojohn_normalization()
df = df["23-47"]

test.display_categories(dataframe=df)

# forestgreen
# salmon
# light-grey
# mediumvioletred
# lightskyblue
# khaki
# cornflowerblue

display1 = test.display_categories(df)
# display2 = cc.apply_style(df=cc.complete_df, columns="23-47")
# display2 = cc.set_for_display(display2)

# display3 = test2.display_categories(df2)
# display4 = cc2.apply_style(df=cc2.complete_df)
# display4 = cc2.set_for_display(display4)

# import pandas as pd
from IPython.display import display_html

# display_html(display3._repr_html_() + display4._repr_html_(), raw=True)
