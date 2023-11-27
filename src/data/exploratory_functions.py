import pandas as pd
import numpy as np
import cobra


class Model_Exploration_Tool:
    def __init__(self, model: cobra.core.model.Model) -> None:
        self.model = model

    def print_reactions_from_metabolite(self, target_metabolite_id: str) -> None:
        target_metabolite = self.model.metabolites.get_by_id(target_metabolite_id)

        print("---- Metabolite Target Information ----")
        print(
            f"ID: {target_metabolite.id}", f"Name: {target_metabolite.name}", sep="\n"
        )
        print("=" * 40)

        for rxn in target_metabolite.reactions:
            print(
                f"ID: {rxn.id}",
                f"Name: {rxn.name}",
                "",
                rxn.build_reaction_string(use_metabolite_names=True),
                rxn.build_reaction_string(use_metabolite_names=False),
                sep="\n",
            )
            print("-" * 40)

    def set_media(
        self,
        medium: list or dict,
        essential=None,
        closed=None,
        lowerbound=0,
        upperbound=1000,
    ):
        for ex_rxn in self.model.exchanges:
            ex_rxn.bounds = (lowerbound, upperbound)

        if type(essential) == list:
            for es_rxn_id in essential:
                es_tgt_rxn = self.model.reactions.get_by_id(es_rxn_id)
                es_tgt_rxn.bounds = (-1000, 1000)

        if type(closed) == list:
            for cl_rxn_id in closed:
                cl_tgt_rxn = self.model.reactions.get_by_id(cl_rxn_id)
                cl_tgt_rxn.bounds = (0, 1000)

        if type(medium) == list:
            for med_comp_rxn_id in medium:
                med_comp_tgt_rxn = self.model.reactions.get_by_id(med_comp_rxn_id)
                med_comp_tgt_rxn.bounds = (-1000, 1000)

        elif type(medium) == dict:
            for med_comp_rxn_id, med_comp_rxn_bounds in medium.items():
                med_comp_tgt_rxn = self.model.reactions.get_by_id(med_comp_rxn_id)
                med_comp_tgt_rxn.bounds = med_comp_rxn_bounds

    #! FVA only works if this the analysis is done on a jupyter notebook file
    def gather_media_fluxes(self, fva_fraction=None):
        up_met_names = {}
        sec_met_names = {}

        # Perform a model FBA optimization otherwise the "summary" function retrieves results of last optimization
        self.model.optimize()
        if fva_fraction != None:
            self.summary_result = self.model.summary(fva_fraction)
        else:
            self.summary_result = self.model.summary()

        self.uptake_df = self.summary_result.uptake_flux.copy()
        self.secretion_df = self.summary_result.secretion_flux.copy()

        for upt_met_id in self.uptake_df["metabolite"]:
            upt_tgt_met = self.model.metabolites.get_by_id(upt_met_id)
            up_met_names[upt_met_id] = upt_tgt_met.name

        self.uptake_df["met. names"] = self.uptake_df["metabolite"].map(up_met_names)
        self.uptake_df = self.uptake_df

        for sec in self.secretion_df["metabolite"]:
            sec_tgt = self.model.metabolites.get_by_id(sec)

            sec_met_names[sec] = sec_tgt.name

        self.secretion_df["met. names"] = self.secretion_df["metabolite"].map(
            sec_met_names
        )
        self.secretion_df = self.secretion_df.loc[
            self.secretion_df["flux"] != 0.00
        ]  # Only interested in those which are predicted to be secreted

        return self.uptake_df, self.secretion_df

    def find_medium_outliers(self, target_df: pd.DataFrame, medium: list or dict):
        if type(medium) == list:
            dif_df = target_df[~target_df["reaction"].isin(medium)]
            return dif_df
        elif type(medium) == dict:
            dif_df = target_df[~target_df["reaction"].isin(medium.keys())]
            return dif_df

    def carbon_source_variation_analysis(self, carbon_list: list) -> pd.DataFrame:
        self.carbon_prediction_results = dict()

        # Optimization without C-Source
        no_c_src_prediction = round(self.model.optimize().objective_value, 5)

        self.carbon_prediction_results["w/o C-Source"] = ["NaN", no_c_src_prediction]

        for c_src_rxn_id in carbon_list:
            with self.model:
                c_src_tgt_rxn = self.model.reactions.get_by_id(c_src_rxn_id)
                c_src_tgt_rxn.bounds = (-1000, 1000)

                c_src_prediction = round(self.model.optimize().objective_value, 6)

                # Summary function must be run before to assign a flux value to the carbon target
                self.model.summary()
                c_src_tgt_flux = round(c_src_tgt_rxn.flux, 5)

                for c_src_met in c_src_tgt_rxn.metabolites:
                    c_src_tgt_met_name = c_src_met.name

                self.carbon_prediction_results[c_src_tgt_rxn.id] = [
                    c_src_tgt_met_name,
                    c_src_prediction,
                    c_src_tgt_flux,
                ]

        self.carbon_prediction_results = pd.DataFrame.from_dict(
            self.carbon_prediction_results,
            orient="index",
            columns=["Met. Name", "FBA Result", "Flux"],
        )

        return self.carbon_prediction_results

    # TODO: Continue from here to correctly name variables and actions within remaining functions
    # TODO: docstring description of each function once finished with optimizing naming

    def reduced_uptake_fba_analysis(self, percentage=1) -> pd.DataFrame:
        # Stores the different pandas series containing the fluxes of the uptaken and secreted metabolites generated in each FBA calculation
        flux_series_dict = dict()
        met_names_dict = dict()
        column_order = ["reaction", "metabolite", "met. names", "flux"]

        # The medium of the model is updated based on which fluxes appear as open (in this case the lower bounds)
        medium = self.model.medium

        # List containing the list of the exchange reaction IDs of the "available" metabolites
        medium_exrxn_list = list(medium.keys())

        self.model.optimize()
        original_solution = self.model.summary().to_frame()
        original_solution = original_solution.loc[
            original_solution["flux"] != 0, ["reaction", "metabolite", "flux"]
        ]

        for met in original_solution["metabolite"]:
            target_met = self.model.metabolites.get_by_id(met)
            met_names_dict[met] = target_met.name

        original_solution["met. names"] = original_solution["metabolite"].map(
            met_names_dict
        )

        original_solution = original_solution.reindex(columns=column_order).rename(
            columns={"flux": "original"}
        )

        # This loop will perform several things:
        # * (Any changes done to the model will not be saved by using the "with" statement)
        ## - Use the available information in the medium to set the reduced uptake values
        ## - Uptake the targets bounds with the reduced value
        ## - Optimize and gather fluxes based on the medium_exrxn_list
        for exrxn in medium_exrxn_list:
            with self.model:
                target_rxn = self.model.reactions.get_by_id(exrxn)
                reduced_uptake_bound = (
                    -medium[exrxn] * percentage
                )  # LB = the uptake bound

                updated_bounds = (reduced_uptake_bound, 1000)
                target_rxn.bounds = updated_bounds

                self.model.optimize()

                solution = self.model.summary().to_frame()
                solution = (
                    solution.loc[solution["flux"] != 0, ["flux"]]
                    .squeeze()
                    .rename(exrxn)
                )

                flux_series_dict[exrxn] = solution

        solution_df = pd.DataFrame(flux_series_dict)

        frames = [original_solution, solution_df]

        complete_solution = pd.concat(frames, axis=1)

        for index in complete_solution.loc[
            complete_solution["reaction"].isnull(), :
        ].index:
            complete_solution.loc[index, "reaction"] = index

            filling_target = self.model.reactions.get_by_id(index)

            for met in filling_target.metabolites:
                complete_solution.loc[index, "metabolite"] = met.id
                complete_solution.loc[index, "met. names"] = met.name

        complete_solution = complete_solution.fillna(0)

        return complete_solution

    def single_varying_uptake_FBA_analysis(self, target_exchange=None):
        # Stores the different pandas series containing the fluxes of the uptaken and secreted metabolites generated in each FBA calculation
        flux_series_dict = dict()
        met_names_dict = dict()
        column_order = ["reaction", "metabolite", "met. names", "flux"]

        ## The medium of the model is updated based on which fluxes appear as open (in this case the lower bounds)
        medium = self.model.medium

        ## List containing the list of the exchange reaction IDs of the "available" metabolites
        medium_exrxn_list = list(medium.keys())

        ## Model needs to be optimized in order to update summary values
        self.model.optimize()

        # Original solution containing base dataframe with values before any change
        original_solution = self.model.summary().to_frame()
        original_solution = original_solution.loc[
            original_solution["flux"] != 0, ["reaction", "metabolite", "flux"]
        ]

        for met in original_solution["metabolite"]:
            target_met = self.model.metabolites.get_by_id(met)
            met_names_dict[met] = target_met.name

        original_solution["met. names"] = original_solution["metabolite"].map(
            met_names_dict
        )

        original_solution = original_solution.reindex(columns=column_order).rename(
            columns={"flux": "original"}
        )

        # This loop will perform several things:
        # * (Any changes done to the model will not be saved by using the "with" statement)
        ## - Use the available information in the medium to set the reduced uptake values
        ## - Uptake the targets bounds with the reduced value
        ## - Optimize and gather fluxes based on the medium_exrxn_list

        if target_exchange != None:
            for exrxn in medium_exrxn_list:
                with self.model:
                    if target_exchange == exrxn:
                        target_rxn = self.model.reactions.get_by_id(exrxn)

                        for percentage in np.arange(0.1, 2.1, 0.1):
                            reduced_uptake_bound = (
                                -medium[exrxn] * percentage
                            )  # LB = the uptake bound
                            updated_bounds = (reduced_uptake_bound, 1000)
                            target_rxn.bounds = updated_bounds

                            self.model.optimize()

                            solution = self.model.summary().to_frame()
                            solution = (
                                solution.loc[solution["flux"] != 0, ["flux"]]
                                .squeeze()
                                .rename(f"FBA {int(percentage*100)}%")
                            )

                            flux_series_dict[f"FBA {int(percentage*100)}%"] = solution

        solution_df = pd.DataFrame(flux_series_dict)

        frames = [original_solution, solution_df]

        complete_solution = pd.concat(frames, axis=1)

        for index in complete_solution.loc[
            complete_solution["reaction"].isnull(), :
        ].index:
            complete_solution.loc[index, "reaction"] = index

            filling_target = self.model.reactions.get_by_id(index)

            for met in filling_target.metabolites:
                complete_solution.loc[index, "metabolite"] = met.id
                complete_solution.loc[index, "met. names"] = met.name

        complete_solution = complete_solution.fillna(0)

        return complete_solution
