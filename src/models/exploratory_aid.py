import cobra
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
from .exploratory_functions import Model_Exploration_Tool as MET
from .path_organizer import PathOrganizer as PO
from pathlib import Path
from datetime import datetime


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
        target_models: list,
        medium: list or dict,
        medium_df: pd.DataFrame,
        closed: list,
        essential: list,
        metabolites: list,
        optimizationType: str = "FBA",
    ):

        LOfModels = self.iter_through_draft_models(
            target_models=target_models,
        )

        biomass_resultsDataframe = pd.DataFrame(index=metabolites["Metabolite"])
        ATPM_resultsDataframe = pd.DataFrame(index=metabolites["Metabolite"])

        for model in LOfModels:

            model: MET

            model.add_and_set_media(
                medium=medium,
                medium_df=medium_df,
                closed_m=closed,
                essential_m=essential,
            )

            biomass_mAssay_results, ATPM_mAssay_results = (
                model.metabolites_growth_energy_optimization_assay(
                    metabolitesDataframe=metabolites,
                    optimizationType=optimizationType,
                )
            )

            biomass_resultsDataframe = biomass_resultsDataframe.join(
                biomass_mAssay_results,
            )
            ATPM_resultsDataframe = ATPM_resultsDataframe.join(
                ATPM_mAssay_results,
            )

        return biomass_resultsDataframe, ATPM_resultsDataframe

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
