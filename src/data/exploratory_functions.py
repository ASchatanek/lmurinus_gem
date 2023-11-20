import pandas as pd


class ExploratoryFunctions:
    def search_reactions(self, model, met_id):
        met_tgt = model.metabolites.get_by_id(met_id)

        for rxn in met_tgt.reactions:
            print(f"{rxn.id} --> {rxn.name}")
            print(rxn.build_reaction_string(use_metabolite_names=True))
            print(rxn.build_reaction_string(use_metabolite_names=False))
            print("")

    def set_media(
        self,
        model,
        medium: list or dict,
        essential=None,
        closed=None,
        lowerbound=0,
        upperbound=1000,
    ):
        for ex in model.exchanges:
            ex.bounds = (lowerbound, upperbound)

        if essential != None and type(essential) == list:
            for es in essential:
                es_tgt = model.reactions.get_by_id(es)
                es_tgt.bounds = (-1000, 1000)

        if closed != None and type(closed) == list:
            for cl in closed:
                x_p = model.reactions.get_by_id(cl)
                x_p.bounds = (0, 1000)

        if type(medium) == list:
            for i in medium:
                i_p = model.reactions.get_by_id(i)
                i_p.bounds = (-1000, 1000)

        elif type(medium) == dict:
            for rxn, bds in medium.items():
                tgt = model.reactions.get_by_id(rxn)
                tgt.bounds = bds

    def gather_media_fluxes(self, model, fva_fraction=None):
        model.optimize()

        if fva_fraction != None:
            self.result = model.summary(fva_fraction)
        else:
            self.result = model.summary()

        uptake_df = self.result.uptake_flux.copy()
        up_met_names = {}

        for upt in uptake_df["metabolite"]:
            upt_tgt = model.metabolites.get_by_id(upt)

            up_met_names[upt] = upt_tgt.name

        uptake_df["met. names"] = uptake_df["metabolite"].map(up_met_names)
        uptake_df = uptake_df

        secretion_df = self.result.secretion_flux.copy()
        sec_met_names = {}

        for sec in secretion_df["metabolite"]:
            sec_tgt = model.metabolites.get_by_id(sec)

            sec_met_names[sec] = sec_tgt.name

        secretion_df["met. names"] = secretion_df["metabolite"].map(sec_met_names)
        secretion_df = secretion_df.loc[
            secretion_df["flux"] != 0.00
        ]  # Only interested in those which were actually secreted

        return uptake_df, secretion_df

    def find_differences(self, target_df, medium: list or dict):
        if type(medium) == list:
            dif_df = target_df[~target_df["reaction"].isin(medium)]
            return dif_df
        elif type(medium) == dict:
            dif_df = target_df[~target_df["reaction"].isin(medium.keys())]
            return dif_df

    def test_carbon_sources(self, model, carbon_list):
        self.carbon_prediction_results = dict()

        # Optimization without C-Source
        carbon_result = round(model.optimize().objective_value, 5)

        self.carbon_prediction_results["w/o C-Source"] = ["NaN", carbon_result]

        for carbon in carbon_list:
            with model:
                carbon_target = model.reactions.get_by_id(carbon)
                carbon_target.bounds = (-1000, 1000)

                carbon_result = round(model.optimize().objective_value, 6)

                # Summary function must be run before to assign a flux value to the carbon target
                model.summary()
                carbon_target_flux = round(carbon_target.flux, 5)

                for metabolite in carbon_target.metabolites:
                    carbon_metabolite_name = metabolite.name

                self.carbon_prediction_results[carbon_target.id] = [
                    carbon_metabolite_name,
                    carbon_result,
                    carbon_target_flux,
                ]

        self.carbon_prediction_results = pd.DataFrame.from_dict(
            self.carbon_prediction_results,
            orient="index",
            columns=["Met. Name", "FBA Result", "Flux"],
        )

        return self.carbon_prediction_results

    def reduced_uptake_fba_analysis(self, model, percentage=1):
        # Stores the different pandas series containing the fluxes of the uptaken and secreted metabolites generated in each FBA calculation
        flux_series_dict = dict()
        met_names_dict = dict()
        column_order = ["reaction", "metabolite", "met. names", "flux"]

        # The medium of the model is updated based on which fluxes appear as open (in this case the lower bounds)
        medium = model.medium

        # List containing the list of the exchange reaction IDs of the "available" metabolites
        medium_exrxn_list = list(medium.keys())

        model.optimize()
        original_solution = model.summary().to_frame()
        original_solution = original_solution.loc[
            original_solution["flux"] != 0, ["reaction", "metabolite", "flux"]
        ]

        for met in original_solution["metabolite"]:
            target_met = model.metabolites.get_by_id(met)
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
            with model:
                target_rxn = model.reactions.get_by_id(exrxn)
                reduced_uptake_bound = (
                    -medium[exrxn] * percentage
                )  # LB = the uptake bound

                updated_bounds = (reduced_uptake_bound, 1000)
                target_rxn.bounds = updated_bounds

                model.optimize()

                solution = model.summary().to_frame()
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

            filling_target = model.reactions.get_by_id(index)

            for met in filling_target.metabolites:
                complete_solution.loc[index, "metabolite"] = met.id
                complete_solution.loc[index, "met. names"] = met.name

        return complete_solution
