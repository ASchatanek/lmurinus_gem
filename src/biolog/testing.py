from scipy.stats import ttest_ind

from notebooks import biolog_variables as bv


bv.abs_590

test1 = bv.abs_590.loc[("23-32", "Plate 1", 168.0)]

test1 = test1.values.tolist()

test2 = bv.abs_590.loc[("23-32", "Plate 2", 168.0)]
test2 = test2.values.tolist()

t_statistic, p_value = ttest_ind(test1, test2)

t_statistic

p_value

import pingouin as pg

result = pg.ttest(test1, test2, correction=True)

result
