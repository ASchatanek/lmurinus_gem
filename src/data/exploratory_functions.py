class ExploratoryFunctions:
    def __init__(self, model):
        self.model = model

    def search_reactions(self, met_id):
        met_tgt = self.model.metabolites.get_by_id(met_id)

        for rxn in met_tgt.reactions:
            print(f"{rxn.id} --> {rxn.name}")
            print(rxn.build_reaction_string(use_metabolite_names=True))
            print("")

    def set_media(
        self,
        medium: list or dict,
        essential: list,
        closed: list,
        lowerbound=0,
        upperbound=1000,
    ):
        for ex in self.model.exchanges:
            ex.bounds = (lowerbound, upperbound)

        if type(medium) == list:
            for i in medium:
                i_p = self.model.reactions.get_by_id(i)
                i_p.bounds = (-1000, 1000)

        elif type(medium) == dict:
            for rxn, bds in medium.items():
                tgt = self.model.reactions.get_by_id(rxn)
                tgt.bounds = bds

        for es in essential:
            es_tgt = self.model.reactions.get_by_id(es)
            es_tgt.bounds = (-1000, 1000)

        for cl in closed:
            x_p = self.model.reactions.get_by_id(cl)
            x_p.bounds = (0, 1000)

    def gather_media_fluxes(self):
        self.model.optimize()
        self.result = self.model.summary()

        uptake_df = self.result.uptake_flux.copy()
        up_met_names = {}

        for upt in uptake_df["metabolite"]:
            upt_tgt = self.model.metabolites.get_by_id(upt)

            up_met_names[upt] = upt_tgt.name

        secretion_df = self.result.secretion_flux.copy()
        sec_met_names = {}

        for sec in secretion_df["metabolite"]:
            sec_tgt = self.model.metabolites.get_by_id(sec)

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
