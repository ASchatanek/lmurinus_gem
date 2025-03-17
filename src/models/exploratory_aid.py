import cobra
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
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
        target_models: list,
        medium: list or dict,
        closed: list,
        essential: list,
        metabolites: list,
    ):

        LOfModels = self.iter_through_draft_models(
            target_models=target_models,
        )

        resultsDf = pd.DataFrame(index=metabolites["Metabolite"])

        for model in LOfModels:

            model: MET

            model.add_missing_medium_met(
                medium=medium,
            )

            model.set_media(
                medium=medium,
                closed=closed,
                essential=essential,
            )

            mAssayResults = model.add_assay_metabolites(
                mets_df=metabolites,
            )

            resultsDf = resultsDf.join(mAssayResults)

        return resultsDf
