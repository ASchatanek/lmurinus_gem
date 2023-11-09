import pandas as pd


class ExploratoryFunctions:
    def search_reactions(self, met_id):
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
        essential: list,
        closed: list,
        lowerbound=0,
        upperbound=1000,
    ):
        for ex in model.exchanges:
            ex.bounds = (lowerbound, upperbound)

        if type(medium) == list:
            for i in medium:
                i_p = model.reactions.get_by_id(i)
                i_p.bounds = (-1000, 1000)

        elif type(medium) == dict:
            for rxn, bds in medium.items():
                tgt = model.reactions.get_by_id(rxn)
                tgt.bounds = bds

        for es in essential:
            es_tgt = model.reactions.get_by_id(es)
            es_tgt.bounds = (-1000, 1000)

        for cl in closed:
            x_p = model.reactions.get_by_id(cl)
            x_p.bounds = (0, 1000)

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

        secretion_df = self.result.secretion_flux.copy()
        sec_met_names = {}

        for sec in secretion_df["metabolite"]:
            sec_tgt = model.metabolites.get_by_id(sec)

            sec_met_names[sec] = sec_tgt.name

        uptake_df["met. names"] = uptake_df["metabolite"].map(up_met_names)
        uptake_df = uptake_df

        secretion_df["met. names"] = secretion_df["metabolite"].map(sec_met_names)
        secretion_df = secretion_df.loc[secretion_df["flux"] != 0.00]

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
