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
        cellStrains = tgt_df.columns.get_level_values(level=0).unique()

        for sample in cellStrains:

            strainTimepoints = (
                tgt_df.loc[:, sample].columns.get_level_values(level=0).unique()
            )

            for tp in strainTimepoints:

                plates = (
                    tgt_df.loc[:, (sample, tp)]
                    .columns.get_level_values(level=0)
                    .unique()
                )

                for plate in plates:
                    tgt = (sample, tp, plate)
                    yield (tgt)

    def determine_average_categories(self, targetData: pd.Series):

        avg_targetData = targetData.mean(axis=1)

        avg_Categories = []

        for avg_Value in avg_targetData:
            if avg_Value < 0.5:
                avg_Categories.append(0)  # 0 is equal to "Indifferent"
            elif avg_Value < 1.5:
                avg_Categories.append(1)  # 1 is equal to "Boundary"
            else:
                avg_Categories.append(2)  # 2 is equal to "Positive"

        return avg_Categories

    def transform_categories_to_intergers(
        self,
        categoriesDataframe: pd.DataFrame,
    ) -> pd.DataFrame:

        def category_to_int(targetCategory):

            # Possible Categories
            indifferent = "I"
            boundary = "B"
            positive = "P"

            # Used to identify any errors while aquiring data
            noValueFound = "E"

            # Assign intergers to possible categories
            if targetCategory == indifferent:
                return 0
            elif targetCategory == boundary:
                return 1
            elif targetCategory == positive:
                return 2
            else:
                return noValueFound

        ## Generate intergerDataframe with categoriesDataframe arquitecture
        intergerDataframe = pd.DataFrame(
            index=categoriesDataframe.index,
            columns=categoriesDataframe.columns,
        )

        ## Iteration
        cellStrains = categoriesDataframe.columns.get_level_values(level=0).unique()
        for strain in cellStrains:
            strainTimepoints = (
                categoriesDataframe.loc[:, strain]
                .columns.get_level_values(level=0)
                .unique()
            )
            for timepoint in strainTimepoints:
                timepointPlates = (
                    categoriesDataframe.loc[:, (strain, timepoint)]
                    .columns.get_level_values(level=0)
                    .unique()
                )

                for plate in timepointPlates:

                    for metabolite in categoriesDataframe.index:

                        targetCombination = (strain, timepoint, plate)

                        metaboliteCategory = categoriesDataframe.loc[
                            metabolite,
                            targetCombination,
                        ]

                        intergerDataframe.loc[metabolite, targetCombination] = (
                            category_to_int(metaboliteCategory)
                        )

        return intergerDataframe

    # * Generate dataframe containing the average category
    def generate_avg_categoriesDataframe(
        self,
        categoriesDataframe: pd.DataFrame,
    ):

        # Generate avg_categoriesDataframe using categoriesDataframe index structure
        avg_CategoriesDataframe = pd.DataFrame(
            index=categoriesDataframe.index,
        )

        # Iteration
        cellStrains = categoriesDataframe.columns.get_level_values(level=0).unique()
        for strain in cellStrains:

            strainTimepoints = (
                categoriesDataframe.loc[:, strain]
                .columns.get_level_values(level=0)
                .unique()
            )

            for timepoint in strainTimepoints:

                tgt_CategoriesDataframe = categoriesDataframe.loc[
                    :, (strain, timepoint)
                ]

                avg_Categories = self.determine_average_categories(
                    targetData=tgt_CategoriesDataframe,
                )

                ## Add resulting category to Results Average Categories df
                avg_CategoriesDataframe[(strain, timepoint)] = avg_Categories

                ## MultiIndex Columns
                ### Generate
                MultiIndexColumns = pd.MultiIndex.from_tuples(
                    avg_CategoriesDataframe.columns,
                    names=["Sample", "TP"],
                )

                ### Assign
                avg_CategoriesDataframe.columns = MultiIndexColumns

                ### Sort
                avg_CategoriesDataframe = avg_CategoriesDataframe.sort_index(axis=1)

        return avg_CategoriesDataframe

    def determine_rangeStatistics(self, targetDataframe: pd.DataFrame):

        # Determine the minimum and maximum value across all plates
        xy_min = targetDataframe.min().min()
        xy_max = targetDataframe.max().max()
        xy_dif = xy_max - xy_min

        # Round up values 3 decimal points
        xy_min = float(round(xy_min, 3))
        xy_max = float(round(xy_max, 3))
        xy_dif = float(round(xy_dif, 3))

        return xy_min, xy_max, xy_dif

    def add_categories_to_targetDataframe(
        self,
        categoriesDataframe: pd.DataFrame,
        targetDataframe: pd.DataFrame,
    ):

        # Generate copy of targetDataframe as merged_targetDataframe
        merged_targetDataframe = targetDataframe.copy()

        cellStrains = categoriesDataframe.columns.get_level_values(level=0).unique()

        for strain in cellStrains:

            strainTimepoints = (
                categoriesDataframe.loc[:, strain]
                .columns.get_level_values(level=0)
                .unique()
            )

            for timepoint in strainTimepoints:

                # Add targetCategories to merged_targetDataframe

                ## Column name of the established categorization
                columnTitle = f"Cat. Mean {str(timepoint)}"

                ## Define targetCategories
                targetCategories = categoriesDataframe.loc[:, (strain, timepoint)]

                ## Add targetCategories
                merged_targetDataframe[(strain, timepoint, columnTitle)] = (
                    targetCategories
                )

                ## Sort merged_targetDataframe
                merged_targetDataframe = merged_targetDataframe.sort_index(axis=1)

        return merged_targetDataframe

    def generate_rangeStatsDataframe(
        self,
        targetDataframe: pd.DataFrame,
    ):

        # Lists
        strainTimepointCombinations = []  # Store (strain, timepoint) combinations
        rangeMinimums = []  # Store minimum value across plates
        rangeMaximums = []  # Story maximum value across plates
        rangeDifferences = []  # Store maximum range

        cellStrains = targetDataframe.columns.get_level_values(level=0).unique()

        for strain in cellStrains:

            strainTimepoints = (
                targetDataframe.loc[:, strain]
                .columns.get_level_values(level=0)
                .unique()
            )

            for timepoint in strainTimepoints:

                # Assemble comb_strainTimepoint
                comb_strainTimepoint = (strain, timepoint)

                # Determine minimumXY, maximumXY and differenceXY
                (
                    min_XY,
                    max_XY,
                    dif_XY,
                ) = self.determine_rangeStatistics(
                    targetDataframe=targetDataframe.loc[:, comb_strainTimepoint]
                )

                # Append information to corresponding lists
                strainTimepointCombinations.append(comb_strainTimepoint)
                rangeMinimums.append(min_XY)
                rangeMaximums.append(max_XY)
                rangeDifferences.append(dif_XY)

        # Generate MultiIndex with comb_strainTimepoints
        comb_strainTimepoint_MultiIndex = pd.MultiIndex.from_tuples(
            strainTimepointCombinations, names=["Sample", "TP"]
        )

        # Establish dictionary for range stats to generate dataframe
        range_dict = {
            "Range Min": rangeMinimums,
            "Range Max": rangeMaximums,
            "Range Dif": rangeDifferences,
        }

        # Range Stats Dataframe
        rangeStatsDataframe = pd.DataFrame(
            range_dict, index=comb_strainTimepoint_MultiIndex
        )

        return rangeStatsDataframe

    def fetch_avg_categoriesDataframe(
        self,
        categoriesDataframe: pd.DataFrame,
    ):

        trans_categoriesDataframe = self.transform_categories_to_intergers(
            categoriesDataframe=categoriesDataframe,
        )

        avg_categoriesDataframe = self.generate_avg_categoriesDataframe(
            categoriesDataframe=trans_categoriesDataframe,
        )

        return avg_categoriesDataframe

    def fetch_merged_targetDataframe(
        self,
        categoriesDataframe,
        targetDataframe,
    ):

        avg_categoriesDataframe = self.fetch_avg_categoriesDataframe(
            categoriesDataframe=categoriesDataframe,
        )

        merged_targetDataframe = self.add_categories_to_targetDataframe(
            categoriesDataframe=avg_categoriesDataframe,
            targetDataframe=targetDataframe,
        )

        return merged_targetDataframe

    def fetch_essentialMetabolites(
        self,
        metabolitesDataframe: pd.DataFrame,
    ):

        essentialMetabolites = [
            [metabolite]
            for metabolite in metabolitesDataframe.loc[
                metabolitesDataframe["Essential"] == 1, "Compound"
            ]
        ]

        return essentialMetabolites

    def generate_plateDistributionsPairGrids(
        self,
        targetData: pd.DataFrame,
        categoriesDataframe: pd.DataFrame,
        strainIDsDataframe: pd.DataFrame,
        kde_thresh: float = 0.1,
        kde_levels: int = 5,
        saveFigures: bool = False,
    ):

        merged_targetDataframe = self.fetch_merged_targetDataframe(
            categoriesDataframe=categoriesDataframe,
            targetDataframe=targetData,
        )

        rangeStatsDataframe = self.generate_rangeStatsDataframe(
            targetDataframe=targetData,
        )

        ##### * INPUTS FOR FIGURES * #####
        ## Input for data type
        dataType = input("What is the target data's type? (BWA, 535nm, 590nm, etc...)")

        ## Input for whether the data has been normalized by BNT
        dataNormStatus = input("Is the target data normalized? (Y, N)")
        if dataNormStatus == "Y":
            dataNormStatus = "normalized"
        else:
            dataNormStatus = "raw"

        ## Input for categorization type
        categoriesType = input("What is the categorization type (BIOLOG, BNT, BEAT)")

        ##### * INPUTS IF saveFigures is TRUE * #####
        ## Input for data type
        if saveFigures is True:
            fileDataName = input(
                "Name of the data dataframe? sample_tp_[input]_categorizationname"
            )

            if fileDataName == "":
                saveFigures = False

        ## Input for categorization type
        if saveFigures is True:
            fileCategoriesName = input(
                "Name of the categorization dataframe? sample_tp_dataname_[input]"
            )

            if fileCategoriesName == "":
                saveFigures = False

        # Iteration
        cellStrains = merged_targetDataframe.columns.get_level_values(level=0).unique()
        for strain in cellStrains:

            # Fetch strain name to substitute ID number
            strainName = strainIDsDataframe.loc[strain, "Strain"]

            strainTimepoints = (
                merged_targetDataframe.loc[:, strain]
                .columns.get_level_values(level=0)
                .unique()
            )
            for timepoint in strainTimepoints:

                # Fetch maximum and minimum from range stats df
                xy_min = rangeStatsDataframe.loc[(strain, timepoint), "Range Min"]
                xy_max = rangeStatsDataframe.loc[(strain, timepoint), "Range Max"]

                # Determine the standard maximum and minimum values to use in graphs
                xy_min = xy_min - (0.2 * xy_max)
                xy_max = xy_max + (0.2 * xy_max)

                # Set current target section of the general mean target dataframe for plotting
                pairGrid_targetData = merged_targetDataframe.loc[:, (strain, timepoint)]

                # Column name of the established categorization
                pairGrid_targetHue = f"Cat. Mean {str(timepoint)}"

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
                    data=pairGrid_targetData,
                    hue=pairGrid_targetHue,
                    corner=False,
                    palette="colorblind",
                    height=4,
                )

                # Name PairGrids by strain name and timepoint
                g.figure.suptitle(
                    f"{strainName} at {timepoint} hrs\nData: {dataType} ({dataNormStatus})   Categorization: {categoriesType} ",
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

                g.map_lower(plot_regplot, color="crimson", data=pairGrid_targetData)

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
                g._legend.set_title("Categories")

                ## Legend Labels
                legendLabels = ["Indifferent", "Boundary", "Positive"]
                for text, label in zip(g._legend.texts, legendLabels):
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
                if saveFigures is True:

                    figureName = f"{str(strain)}_{str(timepoint)}_{str(fileDataName)}_{str(fileCategoriesName)}.png"

                    today_dataFolder = self.generate_daily_folder()
                    figurePath = today_dataFolder / figureName
                    figurePath.resolve()

                    g.savefig(figurePath)

                plt.show()

    # Generate a dataframe which defines range stats and metabolites of interest (positive, boundary) based on average categorization per strain/tp combination
    def generate_summaryPieChartReports(
        self,
        targetData: pd.DataFrame,
        categoriesDataframe: pd.DataFrame,
        strainIDsDataframe: pd.DataFrame,
        metabolitesDataframe: pd.DataFrame,
        saveFigures=False,
    ):

        # Compare lists function
        def highlight_eMets_in_targetList(eMets, targetList):
            highlights = []
            for cval in targetList:
                if cval in eMets:
                    highlights.append(["tab:green"])
                else:
                    highlights.append(["w"])
            return highlights

        # * Series of inputs to name files if saved
        ## Input for data type
        if saveFigures is True:
            fileDataName = input(
                "Name of the data dataframe? sample_tp_[input]_categorizationname"
            )

            if fileDataName == "":
                saveFigures = False

        ## Input for categorization type
        if saveFigures is True:
            fileCategoriesName = input(
                "Name of the categorization dataframe? sample_tp_dataname_[input]"
            )

            if fileCategoriesName == "":
                saveFigures = False

        # * Series of inputs for figures

        ## Input for data type
        dataType = input("What is the data's type? (BWA, 535nm, 590nm, etc...)")

        ## Input for categorization type
        categoriesType = input("What is the categorization type (BIOLOG, BNT, BEAT)")

        merged_targetDataframe = self.fetch_merged_targetDataframe(
            categoriesDataframe=categoriesDataframe,
            targetDataframe=targetData,
        )

        essentialMetabolites = self.fetch_essentialMetabolites(
            metabolitesDataframe=metabolitesDataframe,
        )

        # Iteration
        cellStrains = merged_targetDataframe.columns.get_level_values(level=0).unique()
        for strain in cellStrains:

            # Fetch strain name to substitute ID number
            strainName = strainIDsDataframe.loc[strain, "Strain"]

            strainTimepoints = (
                merged_targetDataframe.loc[:, strain]
                .columns.get_level_values(level=0)
                .unique()
            )
            for timepoint in strainTimepoints:

                avg_categoryColumnTitle = f"Cat. Mean {str(timepoint)}"

                positive_metabolites = [
                    [met]
                    for met in merged_targetDataframe.index
                    if merged_targetDataframe.loc[
                        met, (strain, timepoint, avg_categoryColumnTitle)
                    ]
                    == 2
                ]
                boundary_metabolites = [
                    [met]
                    for met in merged_targetDataframe.index
                    if merged_targetDataframe.loc[
                        met, (strain, timepoint, avg_categoryColumnTitle)
                    ]
                    == 1
                ]
                indifferent_metabolites = [
                    [met]
                    for met in merged_targetDataframe.index
                    if merged_targetDataframe.loc[
                        met, (strain, timepoint, avg_categoryColumnTitle)
                    ]
                    == 0
                ]

                pie_dict = {
                    "Sample - TP": (strain, timepoint),
                    "Strain": strainName,
                    "Positive": positive_metabolites,
                    "Boundary": boundary_metabolites,
                    "Indifferent": indifferent_metabolites,
                }

                ## Creating Datasets
                categoriesLegend = ["Positive", "Boundary", "Indifferent"]
                pieChart_valueCounts = [
                    len(positive_metabolites),
                    len(boundary_metabolites),
                    len(indifferent_metabolites),
                ]

                ## Creating explode pieChart_valueCounts
                pieChart_explode = (0.1, 0.1, 0.1)

                ## Creating color parameters
                pieChart_colors = ("tab:blue", "tab:orange", "tab:gray")

                ## Wedge properties
                pieChart_wedge = {
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
                    pieChart_valueCounts,
                    autopct=lambda pct: autocpt(pct, pieChart_valueCounts),
                    explode=pieChart_explode,
                    shadow=False,
                    colors=pieChart_colors,
                    startangle=90,
                    wedgeprops=pieChart_wedge,
                    textprops=dict(color="black"),
                )

                ## Adding legend
                ax_1.legend(
                    wedges,
                    categoriesLegend,
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
                positive_highlights = highlight_eMets_in_targetList(
                    eMets=essentialMetabolites,
                    targetList=positive_metabolites,
                )
                boundary_hightlights = highlight_eMets_in_targetList(
                    eMets=essentialMetabolites,
                    targetList=boundary_metabolites,
                )

                # Table essential metabolites
                if len(essentialMetabolites) != 0:
                    table_essentialMetabolites = fig.add_subplot(1, 4, 2)
                    table_essentialMetabolites.table(
                        cellText=essentialMetabolites,
                        colLabels=["Essential"],
                        colColours=["tab:green"],
                        loc="upper center",
                        cellLoc="center",
                        edges="closed",
                    ).scale(0.9, 1.3)
                    table_essentialMetabolites.axis("off")

                # Table positive metabolites
                if len(positive_metabolites) != 0:
                    table_positiveMetabolites = fig.add_subplot(1, 4, 3)
                    table_positiveMetabolites.table(
                        cellText=positive_metabolites,
                        colLabels=["Positive"],
                        colColours=["tab:blue"],
                        cellColours=positive_highlights,
                        loc="upper center",
                        cellLoc="center",
                        edges="closed",
                    ).scale(0.9, 1.3)
                    table_positiveMetabolites.axis("off")

                # Table boundary metabolites
                if len(boundary_metabolites) != 0:
                    table_boundaryMetabolites = fig.add_subplot(1, 4, 4)
                    table_boundaryMetabolites.table(
                        cellText=boundary_metabolites,
                        colLabels=["Boundary"],
                        colColours=["tab:orange"],
                        cellColours=boundary_hightlights,
                        loc="upper center",
                        cellLoc="center",
                        edges="closed",
                    ).scale(0.9, 1.3)
                    table_boundaryMetabolites.axis("off")

                fig.suptitle(
                    f"{strainName} at {timepoint} hrs\nData: {dataType}   Categorization: {categoriesType}",
                    y=0.97,
                    fontweight="bold",
                )

                # Save figure
                if saveFigures is True:

                    fig_name = f"sr_{str(strain)}_{str(timepoint)}_{str(fileDataName)}_{str(fileCategoriesName)}.png"

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
        dateString = today.strftime("%Y_%m_%d")

        # Define and resolve the path to the new data folder
        datePath = Path.cwd() / "reports" / "data" / dateString
        datePath.resolve()

        # Generate folder
        Path(datePath).mkdir(
            parents=True,
            exist_ok=True,
        )

        # Return folder's path as a Path object
        return datePath

    ########### GENERATE Summary Categories and Metabolic Assay Table ###########
    def generate_SummaryCategoriesMetAssayTable(
        self,
        categoriesDataframe_BIOLOG: pd.DataFrame,
        categoriesDataframe_BWA: pd.DataFrame,
        categoriesDataframe_590nm: pd.DataFrame,
        categoriesDataframe_750nm: pd.DataFrame,
        strainInfoDataframe: pd.DataFrame,
        biomass_metAssayData: pd.DataFrame,
        mainAPC_metAssayData: pd.DataFrame,
        figureWidth: float = 8,
        figureHeight: float = 20,
        saveFigures: bool = False,
    ):

        # Generate Average Categorization for BIOLOG Categorizations
        avg_CategoriesDataframe_BIOLOG = self.fetch_avg_categoriesDataframe(
            categoriesDataframe=categoriesDataframe_BIOLOG,
        )

        # Generate Average Categorization for BWA Normalized Categorizations
        avg_CategoriesDataframe_BWA = self.fetch_avg_categoriesDataframe(
            categoriesDataframe=categoriesDataframe_BWA,
        )

        # Generate Average Categorization for 590nm Normalized Categorizations
        avg_CategoriesDataframe_590nm = self.fetch_avg_categoriesDataframe(
            categoriesDataframe=categoriesDataframe_590nm,
        )

        avg_CategoriesDataframe_750nm = self.fetch_avg_categoriesDataframe(
            categoriesDataframe=categoriesDataframe_750nm,
        )

        def assign_colors_to_categories(categoriesList):
            categoriesColors = []
            for category in categoriesList:
                if category == 2:
                    categoriesColors.append("tab:blue")
                elif category == 1:
                    categoriesColors.append("tab:orange")
                else:
                    categoriesColors.append("w")

            return categoriesColors

        def assign_colors_to_metAssayData_Biomass(metabolicAssayData):
            metAssayDataColors = []
            for category in metabolicAssayData:
                if category >= 0.1:
                    metAssayDataColors.append("tab:blue")
                elif category < 0.1 and category >= 0.01:
                    metAssayDataColors.append("tab:orange")
                else:
                    metAssayDataColors.append("w")

            return metAssayDataColors

        def assign_colors_to_metAssayData_maintenanceAPC(metAssayData_APC):
            metAssayDataColors_APC = []
            for datapoint in metAssayData_APC:
                if datapoint >= 200:
                    metAssayDataColors_APC.append("tab:blue")
                elif datapoint < 200 and datapoint >= 0.01:
                    metAssayDataColors_APC.append("tab:orange")
                else:
                    metAssayDataColors_APC.append("w")

            return metAssayDataColors_APC

        def arrange_colors_indexwise(
            colorsBIOLOG,
            colorsBWA,
            colors590nm,
            colors750nm,
            colorsBiomassMetAssayData,
            colorsMainAPCMetAssayData,
        ):

            n = len(colorsBIOLOG)
            colorsLists_indexwise = []

            for i in range(n):
                colorsList = []

                colorsList.append(colorsBIOLOG[i])
                colorsList.append(colorsBWA[i])
                colorsList.append(colors590nm[i])
                colorsList.append(colorsMainAPCMetAssayData[i])
                colorsList.append(colors750nm[i])
                colorsList.append(colorsBiomassMetAssayData[i])

                colorsLists_indexwise.append(colorsList)

            return colorsLists_indexwise

        def assign_and_arrange_colors_to_summaryDataframe():

            colorsBIOLOG = assign_colors_to_categories(summaryDataframe["BIOLOG"])
            colorsBWA = assign_colors_to_categories(summaryDataframe["BWA"])
            colors590nm = assign_colors_to_categories(summaryDataframe["590nm"])
            colors750nm = assign_colors_to_categories(summaryDataframe["750nm"])
            colorsBiomassMetAssayData = assign_colors_to_metAssayData_Biomass(
                summaryDataframe["BioM"]
            )
            colorsMainAPCMetAssayData = assign_colors_to_metAssayData_maintenanceAPC(
                summaryDataframe["mainAPC"]
            )

            colors = arrange_colors_indexwise(
                colorsBIOLOG=colorsBIOLOG,
                colorsBWA=colorsBWA,
                colors590nm=colors590nm,
                colors750nm=colors750nm,
                colorsBiomassMetAssayData=colorsBiomassMetAssayData,
                colorsMainAPCMetAssayData=colorsMainAPCMetAssayData,
            )

            return colors

        # Iteration
        cellStrains = avg_CategoriesDataframe_BIOLOG.columns.get_level_values(
            level=0
        ).unique()
        for strain in cellStrains:

            # Fetch strain name to substitute ID number
            strainName = strainInfoDataframe.loc[strain, "Strain"]
            strainDSMZ = strainInfoDataframe.loc[strain, "DSMZ-number"]

            targetBiomass_metAssayData = biomass_metAssayData.loc[:, strainDSMZ]
            targetMainAPC_metAssayData = mainAPC_metAssayData.loc[:, strainDSMZ]

            targetBiomass_metAssayData = targetBiomass_metAssayData.round(2)
            targetMainAPC_metAssayData = targetMainAPC_metAssayData.round(2)

            strainTimepoints = (
                avg_CategoriesDataframe_BIOLOG.loc[:, strain]
                .columns.get_level_values(level=0)
                .unique()
            )
            for timepoint in strainTimepoints:

                targetCategories_BIOLOG = avg_CategoriesDataframe_BIOLOG.loc[
                    :, (strain, timepoint)
                ]
                targetCategories_BWA = avg_CategoriesDataframe_BWA.loc[
                    :, (strain, timepoint)
                ]
                targetCategories_590nm = avg_CategoriesDataframe_590nm.loc[
                    :, (strain, timepoint)
                ]
                targetCategories_750nm = avg_CategoriesDataframe_750nm.loc[
                    :, (strain, timepoint)
                ]

                summaryDataframe = pd.DataFrame(
                    index=avg_CategoriesDataframe_BIOLOG.index,
                )

                summaryDataframe[("BIOLOG")] = targetCategories_BIOLOG
                summaryDataframe[("BWA")] = targetCategories_BWA
                summaryDataframe[("590nm")] = targetCategories_590nm
                summaryDataframe[("mainAPC")] = targetMainAPC_metAssayData
                summaryDataframe[("750nm")] = targetCategories_750nm
                summaryDataframe[("BioM")] = targetBiomass_metAssayData

                summaryDataframe = summaryDataframe.sort_values(
                    summaryDataframe.columns.tolist(),
                    ascending=False,
                )

                summaryColors = assign_and_arrange_colors_to_summaryDataframe()

                # Table Generation with Categories and Metabolic Assay Data and Colors
                fig, ax = plt.subplots(figsize=(figureWidth, figureHeight))

                summaryCategoriesMetAssayTable = ax.table(
                    cellText=summaryDataframe.values,
                    cellColours=summaryColors,
                    rowLabels=summaryDataframe.index,
                    rowLoc="right",
                    colLabels=summaryDataframe.columns,
                    cellLoc="center",
                    loc="best",
                )

                summaryCategoriesMetAssayTable.auto_set_font_size(False)
                summaryCategoriesMetAssayTable.set_fontsize(12)

                fig.suptitle(f"{strainName} : {timepoint} hrs")
                ax.axis("off")
                fig.tight_layout()

                # Save figure
                if saveFigures is True:

                    figureFileName = f"{str(strain)}_{str(timepoint)}_sumtable.png"

                    todaysDateFolder = self.generate_daily_folder()
                    figureFilePath = todaysDateFolder / figureFileName
                    figureFilePath.resolve()

                    fig.savefig(figureFilePath)

                plt.show()
