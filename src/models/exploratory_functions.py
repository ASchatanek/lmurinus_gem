import cobra
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class Model_Exploration_Tool:

    """Exploration tool for GEM models. Includes a variety of methods to display reactions, set media components and analyse media contraining effects on the model´s predictions."""

    def __init__(self, model: cobra.core.model.Model) -> None:
        self.model = model
        self.model.optimize()

    def print_reactions_from_metabolite(self, target_metabolite_id: str) -> None:
        """Prints reactions associated to a specific metabolite ID. It displayes the individual reactions´ ID, name along with the reactions stochiometry with IDs and names.

        Parameters
        ----------
        target_metabolite_id : str
            Target metabolite ID for which the reactions should be printed out.
        """
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
        essential: list = None,
        closed: list = None,
        lowerbound: float = 0.0,
        upperbound: float = 1000.0,
    ):
        """Updates the current model´s media by fetching exchange reaction IDs and their bounds from different dictionaries or lists.

        Parameters
        ----------
        medium : listordict
            Exchange reaction IDs for the media components stored as a dictionary or list. In dictionaries the IDs and their bounds (as tuple) should be stored as key and values respectively. Bounds for IDs contained in lists are set as (-1000,1000).
        essential : list, optional
            List of exchange reaction IDs that are set as open, meaning their bounds are (-1000, 1000). Default is None.
        closed : list, optional
            List of exchange reaction IDs that are set as closed for uptake, meaning their bounds are (0, 1000). Default is None.
        lowerbound : float, optional
            Float value that define the lowerbound for any exchange reaction that isn´t included in "medium", "essential" and "closed" variables. Default is 0.0.
        upperbound : float, optional
            Float value that define the upperbound for any exchange reaction that isn´t included in "medium", "essential" and "closed" variables. Default is 1000.0.
        """
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
    def gather_media_fluxes(self, fva_fraction: float = None) -> pd.DataFrame:
        """Performs a FBA optimization on the model and returns dataframes containing the resulting uptaken and secreted metabolites along with their flux information.

        Parameters
        ----------
        fva_fraction : float, optional
            Decimal fraction representing the FVA optimization percentage. If defined the dataframes will also contain resulting FVA predictions. WARNING: This attribute only works if the method is done on a jupyter notebook.

        Returns
        -------
        pd.Dataframe
            Returns the uptaken and secreted metabolites with their predicted fluxes (and FVA predictions if defined) as a pandas Dataframe.
        """
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

    def find_medium_outliers(
        self, target_df: pd.DataFrame, medium: list or dict
    ) -> pd.DataFrame:
        """Identifies any metabolites that result from an FBA optimization and are not included in the current model´s medium.

        Parameters
        ----------
        target_df : pd.DataFrame
            Pandas dataframe in from which possible outliers should be identified from. Intented to be used with the resulting dataframes from method "gather_media_fluxes".
        medium : listordict
            List of dictionary containing the exchange reaction IDs. To be used to identify any exchange reaction not included in this list but included in the target_df variable. Intended to be used with the medium list of dictionary used in method "set_media".

        Returns
        -------
        pd.DataFrame
            Returns dataframe containing any outliers and their information not present in the medium list/dictionary but included in the target_df dataframe.
        """
        if type(medium) == list:
            dif_df = target_df[~target_df["reaction"].isin(medium)]
            return dif_df
        elif type(medium) == dict:
            dif_df = target_df[~target_df["reaction"].isin(medium.keys())]
            return dif_df

    def carbon_source_variation_analysis(self, carbon_list: list) -> pd.DataFrame:
        """Generates a pandas dataframe containing the resulting FBA predictions for different carbon sources on the current model´s media as well as the prediction without any source.

        Parameters
        ----------
        carbon_list : list
            List containing the exchange reaction IDs for different carbon sources.

        Returns
        -------
        pd.DataFrame
            Returns a pandas dataframe containing the carbon source metabolite name, the predicted objective function value and the predicted flux for their corresponding exchange reaction.
        """
        c_source_prediction_results = dict()

        # Optimization without C-Source
        no_c_src_prediction = round(self.model.optimize().objective_value, 5)

        c_source_prediction_results["w/o C-Source"] = ["NaN", no_c_src_prediction]

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

                c_source_prediction_results[c_src_tgt_rxn.id] = [
                    c_src_tgt_met_name,
                    c_src_prediction,
                    c_src_tgt_flux,
                ]

        self.c_source_prediction_results_df = pd.DataFrame.from_dict(
            c_source_prediction_results,
            orient="index",
            columns=["Met. Name", "FBA Result", "Flux"],
        )

        return self.c_source_prediction_results_df

    # TODO: Continue from here to correctly name variables and actions within remaining functions
    # TODO: docstring description of each function once finished with optimizing naming

    def reduced_uptake_fba_analysis(self, percentage: float = 1.0) -> pd.DataFrame:
        """Method to analyse the effects of contraining the uptake of individual medium components on the model´s predicted secretion and uptake fluxes. An iterative process takes the individual components of the current defined model´s medium. Then, while using the "with" statement, modifies the current uptake bound (neg. lowerbound) by the defined percentage and performs an FBA optimization. Using the "summary()" method, the secreted and uptaked fluxes are retrieved and stored.

        Parameters
        ----------
        percentage : float, optional
            Decimal value (1.0 being 100%) to which the individual components uptake bound should be modified to. Default is 1.0.

        Returns
        -------
        pd.DataFrame
            Returns a pandas dataframe containing the secretion and uptake resulting fluxes for each contrained component of the current model´s medium. Additionally, an "original" column containing the predicted secretion and uptake results without any contrained uptake bounds is added as reference.
        """

        # Stores the different pandas series containing the fluxes of the uptaken and secreted metabolites generated in each FBA calculation
        flux_series_dict = dict()
        met_names_dict = dict()
        column_order = ["reaction", "metabolite", "met. names", "flux"]

        # The medium of the model is updated based on which fluxes appear as open (in this case the lower bounds)
        medium = self.model.medium

        # List containing the exchange reaction IDs of the "available" metabolites
        medium_ex_rxn_list = list(medium.keys())

        self.model.optimize()
        original_solution = self.model.summary().to_frame()
        original_solution = original_solution.loc[
            original_solution["flux"] != 0, ["reaction", "metabolite", "flux"]
        ]

        for met_id in original_solution["metabolite"]:
            tgt_met = self.model.metabolites.get_by_id(met_id)
            met_names_dict[met_id] = tgt_met.name

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
        ## - Optimize and gather fluxes based on the medium_ex_rxn_list
        for ex_rxn_id in medium_ex_rxn_list:
            with self.model:
                tgt_rxn = self.model.reactions.get_by_id(ex_rxn_id)
                reduced_uptake_bound = (
                    -medium[ex_rxn_id] * percentage
                )  # LB = the uptake bound

                updated_bounds = (reduced_uptake_bound, 1000)
                tgt_rxn.bounds = updated_bounds

                self.model.optimize()

                solution = self.model.summary().to_frame()
                solution = (
                    solution.loc[solution["flux"] != 0, ["flux"]]
                    .squeeze()
                    .rename(ex_rxn_id)
                )

                flux_series_dict[ex_rxn_id] = solution

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

    def fetch_constrained_medium_fba_fluxes(
        self, percentage: float = 1.0
    ) -> pd.DataFrame:
        """Gathers resulting fluxes of the current medium components after perfoming an optimization while iteratively constraining one individual components´ bounds to their predicted optimum flux.

        Parameters
        ----------
        percentage : float, optional
            Percentile decimal value to which the predicted optimum flux of each medium component should be contrained to. Default is 1.0 (100%).

        Returns
        -------
        pd.DataFrame
            Returns dataframe containing the predicted fluxes for the components of the model´s current medium.
        """
        # List containing the list of the exchange reaction IDs of the "available" metabolites
        column_order = ["reaction", "metabolite", "met. names", "flux"]
        medium_ex_rxn_list = list(self.model.medium.keys())
        met_info_dict = dict()
        flux_series_dict = dict()

        self.model.optimize()
        reference_solution = (
            self.model.summary()
            .to_frame()
            .loc[medium_ex_rxn_list, ["reaction", "metabolite", "flux"]]
        )

        for met_id in reference_solution["metabolite"]:
            tgt_met = self.model.metabolites.get_by_id(met_id)
            met_info_dict[met_id] = tgt_met.name

        reference_solution["met. names"] = reference_solution["metabolite"].map(
            met_info_dict
        )

        reference_solution = reference_solution.reindex(columns=column_order).rename(
            columns={"flux": "original"}
        )

        reference_solution["exchange type"] = reference_solution["original"].apply(
            self.__define_flux_type
        )

        for ex_rxn_id in medium_ex_rxn_list:
            with self.model:
                tgt_rxn = self.model.reactions.get_by_id(ex_rxn_id)

                tgt_rxn_flux = reference_solution.loc[ex_rxn_id, "original"]

                reduced_contrained_bound = -tgt_rxn_flux * percentage

                updated_bounds = (reduced_contrained_bound, reduced_contrained_bound)

                tgt_rxn.bounds = updated_bounds

                self.model.optimize()

                solution = self.model.summary().to_frame()
                solution = (
                    solution.loc[medium_ex_rxn_list, ["flux"]]
                    .squeeze()
                    .rename(ex_rxn_id)
                )

                flux_series_dict[ex_rxn_id] = solution

        solution_df = pd.DataFrame(flux_series_dict)

        frames = [reference_solution, solution_df]

        complete_solution = pd.concat(frames, axis=1)

        return complete_solution

    def fetch_constrained_sec_and_upt_fluxes(
        self, percentage: float = 1.0
    ) -> pd.DataFrame:
        # List containing the list of the exchange reaction IDs of the "available" metabolites
        column_order = ["reaction", "metabolite", "met. names", "flux"]
        medium_ex_rxn_list = list(self.model.medium.keys())
        met_info_dict = dict()
        flux_series_dict = dict()

        self.model.optimize()
        reference_solution = self.model.summary().to_frame()
        reference_solution = reference_solution.loc[
            reference_solution["flux"] != 0, ["reaction", "metabolite", "flux"]
        ]

        for met_id in reference_solution["metabolite"]:
            tgt_met = self.model.metabolites.get_by_id(met_id)
            met_info_dict[met_id] = tgt_met.name

        reference_solution["met. names"] = reference_solution["metabolite"].map(
            met_info_dict
        )

        reference_solution = reference_solution.reindex(columns=column_order).rename(
            columns={"flux": "original"}
        )

        reference_solution["exchange type"] = reference_solution["original"].apply(
            self.__define_flux_type
        )

        for ex_rxn_id in medium_ex_rxn_list:
            if ex_rxn_id in reference_solution.index:
                with self.model:
                    tgt_rxn = self.model.reactions.get_by_id(ex_rxn_id)

                    tgt_rxn_flux = reference_solution.loc[ex_rxn_id, "original"]

                    reduced_contrained_bound = -tgt_rxn_flux * percentage

                    updated_bounds = (
                        reduced_contrained_bound,
                        reduced_contrained_bound,
                    )

                    tgt_rxn.bounds = updated_bounds

                    self.model.optimize()

                    solution = self.model.summary().to_frame()
                    solution = (
                        solution.loc[reference_solution.index, ["flux"]]
                        .squeeze()
                        .rename(ex_rxn_id)
                    )

                    flux_series_dict[ex_rxn_id] = solution

        solution_df = pd.DataFrame(flux_series_dict)

        frames = [reference_solution, solution_df]

        complete_solution = pd.concat(frames, axis=1)

        return complete_solution

    def __define_flux_type(self, value):
        if value < 0:
            return "Secretion"
        if value == 0:
            return "No Flux"
        elif value > 0:
            return "Uptake"

    def constrained_fba_analysis(self, percentage: float = 1.0) -> pd.DataFrame:
        complete_constr_df = self.fetch_constrained_sec_and_upt_fluxes(
            percentage=percentage
        )

        sec_complete_df = complete_constr_df[
            complete_constr_df["exchange type"] == "Secretion"
        ]
        sec_results_df = sec_complete_df.iloc[:, 5:].copy()

        upt_complete_df = complete_constr_df[
            complete_constr_df["exchange type"] == "Uptake"
        ]
        upt_results_df = upt_complete_df.iloc[:, 5:].copy()

        for column in sec_results_df:
            sec_results_df[column] = (
                sec_results_df[column] / sec_complete_df["original"]
            )
            sec_results_df[column] = sec_results_df[column]

        sec_results_df = sec_results_df.round(2)

        for column in upt_results_df:
            upt_results_df[column] = (
                upt_results_df[column] / upt_complete_df["original"]
            )

        upt_results_df = upt_results_df.round(2)

        return sec_results_df, upt_results_df

    def display_hm_constrained_medium_fba_analysis(self, percentage: float = 1.0):
        sec_df, upt_df = self.constrained_fba_analysis(percentage=percentage)

        fig, ax = plt.subplots(figsize=(25, 20))
        outliers = upt_df.map(lambda v: v if v > 1 or v < -1 else "")
        upt_hm = sns.heatmap(
            data=upt_df,
            linewidths=0.1,
            ax=ax,
            annot=outliers,
            annot_kws={"fontsize": 9},
            fmt="",
            cmap="PRGn",
            vmax=2,
            vmin=-2,
        )

        upt_x_labels = []
        upt_y_labels = []

        for column in upt_hm.get_xticklabels():
            c_text = column.get_text()
            target_rxn = self.model.reactions.get_by_id(c_text)
            for met in target_rxn.metabolites:
                target_met = met.name
                upt_x_labels.append(target_met)

        for index in upt_hm.get_yticklabels():
            i_text = index.get_text()
            target_rxn = self.model.reactions.get_by_id(i_text)
            for met in target_rxn.metabolites:
                target_met = met.name
                upt_y_labels.append(target_met)

        upt_hm.set_xticklabels(upt_x_labels, rotation=45, horizontalalignment="left")
        upt_hm.set_yticklabels(upt_y_labels)
        upt_hm.xaxis.tick_top()
        plt.tight_layout()
        plt.show()

        fig, ax = plt.subplots(figsize=(25, 20))
        outliers = sec_df.map(lambda v: v if v > 1 or v < -1 else "")
        sec_hm = sns.heatmap(
            data=sec_df,
            linewidths=0.1,
            ax=ax,
            annot=outliers,
            annot_kws={"fontsize": 9},
            fmt="",
            cmap="PRGn",
            vmax=2,
            vmin=-2,
        )

        sec_x_labels = []
        sec_y_labels = []

        for column in sec_hm.get_xticklabels():
            c_text = column.get_text()
            target_rxn = self.model.reactions.get_by_id(c_text)
            for met in target_rxn.metabolites:
                target_met = met.name
                sec_x_labels.append(target_met)

        for index in sec_hm.get_yticklabels():
            i_text = index.get_text()
            target_rxn = self.model.reactions.get_by_id(i_text)
            for met in target_rxn.metabolites:
                target_met = met.name
                sec_y_labels.append(target_met)

        sec_hm.set_xticklabels(sec_x_labels, rotation=45, horizontalalignment="left")
        sec_hm.set_yticklabels(sec_y_labels)
        sec_hm.xaxis.tick_top()
        plt.tight_layout()
        plt.show()

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
