import re
from datetime import datetime
from pathlib import Path

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .exploratory_functions import Model_Exploration_Tool as MET
from .path_organizer import PathOrganizer as PO


class ExploratoryAid:

    def __init__(self):
        pass

    def iter_through_draft_models(self, target_models: list):

        list_of_models = []

        for model in target_models:

            # Path Organizer
            po = PO()

            # target draft model directory
            tgtmDir = po.get_draft_model_path(draft_name=model)

            # load target draft model
            tgtm = po.load_model(tgtmDir)

            # load MET for target draft model
            tgtMET = MET(model=tgtm)

            # attach target model
            list_of_models.append(tgtMET)

        return list_of_models

    def metabolite_assay(
        self,
        targetModels: list,
        mediumMetabolites: list or dict,
        mediumMetabolitesDataframe: pd.DataFrame,
        closedMetabolites: list,
        targetModels_EssentialMetabolites: dict,
        assayMetabolitesDataframe: pd.DataFrame,
        optimizationType: str = "FBA",
    ):

        LOfModels = self.iter_through_draft_models(
            target_models=targetModels,
        )

        biomass_resultsDataframe = pd.DataFrame(
            index=assayMetabolitesDataframe["Metabolite"]
        )
        ATPM_resultsDataframe = pd.DataFrame(
            index=assayMetabolitesDataframe["Metabolite"]
        )

        metaboliteExistence_resultsDataframe = pd.DataFrame(
            index=assayMetabolitesDataframe["Metabolite"]
        )

        for model in LOfModels:

            model: MET

            target_essentialMetabolites = targetModels_EssentialMetabolites[
                model.model.id
            ]

            model.add_and_set_media(
                medium=mediumMetabolites,
                medium_df=mediumMetabolitesDataframe,
                closed_m=closedMetabolites,
                essential_m=target_essentialMetabolites,
            )

            (
                biomass_mAssay_results,
                ATPM_mAssay_results,
                MetabolitesExistence_results,
            ) = model.metabolites_growth_energy_optimization_assay(
                metabolitesDataframe=assayMetabolitesDataframe,
                optimizationType=optimizationType,
            )

            biomass_resultsDataframe = biomass_resultsDataframe.join(
                biomass_mAssay_results,
            )

            ATPM_resultsDataframe = ATPM_resultsDataframe.join(
                ATPM_mAssay_results,
            )

            metaboliteExistence_resultsDataframe = (
                metaboliteExistence_resultsDataframe.join(
                    MetabolitesExistence_results,
                )
            )

        return (
            biomass_resultsDataframe,
            ATPM_resultsDataframe,
            metaboliteExistence_resultsDataframe,
        )

    def metabolite_assay_testing(
        self,
        targetModels: list,
        mediumMetabolites: list or dict,
        mediumMetabolitesDataframe: pd.DataFrame,
        closedMetabolites: list,
        targetModels_EssentialMetabolites: dict,
        assayMetabolitesDataframe: pd.DataFrame,
        optimizationType: str = "FBA",
    ):

        results_metAssays = pd.DataFrame()

        LOfModels = self.iter_through_draft_models(
            target_models=targetModels,
        )

        for model in LOfModels:

            model: MET

            target_essentialMetabolites = targetModels_EssentialMetabolites[
                model.model.id
            ]

            model.add_and_set_media(
                medium=mediumMetabolites,
                medium_df=mediumMetabolitesDataframe,
                closed_m=closedMetabolites,
                essential_m=target_essentialMetabolites,
            )

            model_results_metAssay = model.testing_something(
                metabolitesDataframe=assayMetabolitesDataframe,
                optimizationType=optimizationType,
            )

            results_metAssays = pd.concat(
                [results_metAssays, model_results_metAssay],
                axis=1,
            )

        return results_metAssays

    def generate_essentialMetabolites_dictionaries(
        self,
        targetModels: list,
        mediumMetabolites: list or dict,
        mediumMetabolitesDataframe: pd.DataFrame,
        closedMetabolites: list,
    ):

        # Define empty result dictionaries
        dict_crucial_eMetabolites = dict()
        dict_variable_eMetabolites = dict()

        LOfModels = self.iter_through_draft_models(
            target_models=targetModels,
        )

        for model in LOfModels:

            model: MET

            crucial_essentialMetabolites, variable_essentialMetabolites = (
                model.find_crucial_and_variable_essentialMetabolites(
                    baseMedium=mediumMetabolites,
                    baseMediumDataframe=mediumMetabolitesDataframe,
                    closedMetabolites=closedMetabolites,
                )
            )

            dict_crucial_eMetabolites[model.model.id] = crucial_essentialMetabolites
            dict_variable_eMetabolites[model.model.id] = variable_essentialMetabolites

        return dict_crucial_eMetabolites, dict_variable_eMetabolites

    def essentialMetabolites_report(
        self,
        targetModels: list,
        mediumMetabolites: list or dict,
        mediumMetabolitesDataframe: pd.DataFrame,
        closedMetabolites: list,
        crucial_eMetabolites: dict,
        variable_eMetabolites: dict,
    ):

        LOfModels = self.iter_through_draft_models(target_models=targetModels)

        best_essentialMetabolites = dict()

        # Iterate through models
        for model in LOfModels:

            model: MET

            # Retrieve target crucial and variable essentialMetabolites
            model_crucial_eMets = crucial_eMetabolites[model.model.id]
            model_variable_eMets = variable_eMetabolites[model.model.id]

            # Print for crucial essential without any extra variable essential metabolite
            model.add_and_set_media(
                medium=mediumMetabolites,
                medium_df=mediumMetabolitesDataframe,
                closed_m=closedMetabolites,
                essential_m=model_crucial_eMets,
            )

            slimResult = model.model.slim_optimize()

            print(f"{model.model.id}")
            print("......................")
            print("Essential Metabolites:")
            print("")
            for crucial in model_crucial_eMets:

                crucial_name = model.model.reactions.get_by_id(crucial).name

                print(f"{crucial} --> {crucial_name}")

            print("")
            print(f"FBA: {round(slimResult, 4)}")
            print("")
            if not model_variable_eMets:

                best_essentialMetabolites[model.model.id] = model_crucial_eMets

            else:

                best_essentialMetabolites_result = 0.0

                print("......................")
                print("Variable Metabolites:")
                print("")
                for variable in model_variable_eMets:

                    essentialMetabolites = model_crucial_eMets + [variable]

                    variable_name = model.model.reactions.get_by_id(variable).name

                    model.set_media(
                        medium=mediumMetabolites,
                        essential=essentialMetabolites,
                        closed=closedMetabolites,
                    )

                    slimResult = model.model.slim_optimize()

                    print(f"{variable} --> {variable_name}")
                    print(f"FBA: {round(slimResult, 4)}")
                    print("")

                    if slimResult > best_essentialMetabolites_result:

                        best_essentialMetabolites_result = slimResult
                        best_essentialMetabolites_combination = essentialMetabolites

                best_essentialMetabolites[model.model.id] = (
                    best_essentialMetabolites_combination
                )

            print("_______________________________")

        return best_essentialMetabolites

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

    def generate_daily_pickle_folder(self) -> Path:

        todaysFolderPath = self.generate_daily_folder()

        todaysPickleFolderPath = todaysFolderPath / "pickles"

        todaysPickleFolderPath.resolve()

        # Generate folder
        Path(todaysPickleFolderPath).mkdir(
            parents=True,
            exist_ok=True,
        )

        return todaysPickleFolderPath

    def fetch_pickleFilePath(self) -> Path:

        todaysPickleFolderPath = self.generate_daily_pickle_folder()

        pickleFileName = input("File name?")

        pickleFilePath = todaysPickleFolderPath / f"{pickleFileName}.pkl"
        pickleFilePath.resolve()

        return pickleFilePath

    def save_dataframe_to_pickle(self, dataframe: pd.DataFrame):

        pickleFilePath = self.fetch_pickleFilePath()

        dataframe.to_pickle(pickleFilePath)

    def load_dataframe_from_pickle(self) -> pd.DataFrame:

        pickleFilePath = self.fetch_pickleFilePath()

        dataframeFromPickle = pd.read_pickle(pickleFilePath)

        return dataframeFromPickle
