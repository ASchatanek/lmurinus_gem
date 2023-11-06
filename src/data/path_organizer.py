from pathlib import Path
from cobra.io import read_sbml_model, write_sbml_model
import os
import re


class PathOrganizer:
    def __init__(self):
        self.data_dir = Path()
        self.model_dir = Path()
        self.model_version_dir = Path()
        self.target_mo_version_dir = Path()
        self.paths = {}
        self.mo_name = None
        self.major_version = None
        self.minor_version = None
        self.sv_pattern = rf"\d+\.\d+\.(\d+)\.sbml$"

    # Functions
    def get_draft_model_path(self, draft_name):
        self.data_dir = Path.cwd() / "data" / "raw" / draft_name
        self.data_dir = self.data_dir.resolve()
        return self.data_dir

    def get_models_main_folder_path(self, mo_name=None):
        if self.mo_name != None and mo_name == None:
            self.model_dir = Path.cwd() / "models" / self.mo_name
        else:
            if mo_name != None:
                self.mo_name = mo_name
                self.model_dir = Path.cwd() / "models" / self.mo_name
            else:
                self.mo_name = str(input("Name of Organism?"))
                self.model_dir = Path.cwd() / "models" / self.mo_name

        self.model_dir = self.model_dir.resolve()
        if self.model_dir.is_dir() is False:
            print("Organism not found.")
            self.mo_name = None

        else:
            return self.model_dir

    def get_model_version_path(self, mo_name=None):
        self.get_models_main_folder_path(mo_name=mo_name)

        self.target_mo_version = None

        print("Available Model Versions:")
        for number, model_id in enumerate(os.listdir(self.model_dir)):
            print(f"{number}) {model_id}")

        self.target_mo_version = str(input("What model should be loaded?"))

        match = re.match(self.sv_pattern, self.target_mo_version)

        if match:
            self.target_mo_version_dir = os.path.join(
                self.model_dir, self.target_mo_version
            )
            return self.target_mo_version_dir
        else:
            "Model not found."

    def load_model(self, path):
        self.model = read_sbml_model(str(path))
        return self.model

    def get_existing_models(self):
        self.existing_models = []

        self.version_pattern = (
            rf"{self.major_version}\.{self.minor_version}\.(\d+)\.sbml$"
        )

        for model_id in os.listdir(self.model_dir):
            match = re.match(self.version_pattern, model_id)
            if match:
                self.existing_models.append(model_id)

        return self.existing_models

    def semantic_versioning(self):
        if self.model_dir.is_dir() is False:
            self.model_dir = self.get_models_main_folder_path()

        self.get_existing_models()

        if self.existing_models:
            self.patch_versions = [
                int(re.match(self.sv_pattern, model_id).group(1))
                for model_id in self.existing_models
            ]
            self.next_patch = max(self.patch_versions) + 1
        else:
            self.next_patch = 1

        self.model_filename = (
            f"{self.major_version}.{self.minor_version}.{self.next_patch}.sbml"
        )

    def save_model(self):
        self.get_models_main_folder_path()

        if self.major_version is None:
            self.major_version = int(input("Major Version Number?"))

        if self.minor_version is None:
            self.minor_version = int(input("Minor Version Number?"))

        self.semantic_versioning()

        print("Existing Model Versions:")
        for number, model_id in enumerate(self.existing_models):
            print(f"{number}) {model_id}")
        print("-" * 50)

        print(f"Model version will be saved as {self.model_filename}")

        sure = input("Are you sure? (y/n)")

        if sure == "y":
            self.model_version_dir = os.path.join(self.model_dir, self.model_filename)
            write_sbml_model(cobra_model=self.model, filename=self.model_version_dir)
            print("-" * 50)
            print("Version saved.")
        else:
            print("-" * 50)
            print("Nothing was saved.")
