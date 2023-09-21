from pathlib import Path
from cobra.io import load_model, read_sbml_model, write_sbml_model


class PathOrganizer:
    def __init__(self):
        self.data_dir = Path()
        self.xml_path = ""
        self.sbml_path = ""
        self.files_path_list = list()
        self.paths = {}

        # Initiate Functions
        # self.get_main_folder_path()
        # self.get_and_organize_files_path()

    # Functions
    def get_draft_model_path(self, draft_name):
        self.data_dir = Path.cwd() / "data" / "raw" / draft_name
        self.data_dir = self.data_dir.resolve()
        return str(self.data_dir)

    def get_models_main_folder_path(self, target_mo_name):
        self.data_dir = Path.cwd() / "models" / target_mo_name
        self.data_dir = self.data_dir.resolve()
        return self.data_dir

    def get_and_organize_files_path(self):
        xml_suf = ".xml"
        sbml_suf = ".sbml"

        xml_files = []
        sbml_files = []

        for file in self.data_dir.glob("*.*"):
            file = file.resolve()
            self.files_path_list.append(file)

        for path in self.files_path_list:
            if str(path).endswith(xml_suf):
                xml_files.append(path)

            if str(path).endswith(sbml_suf):
                sbml_files.append(path)

        self.paths["XML"] = xml_files
        self.paths["SBML"] = sbml_files

    def define_xml_path(self):
        options = self.paths["XML"]
        for i, item in enumerate(options):
            print(i, "||", item)
        x = int(input("Which model do you want? Start with 0: "))
        self.xml_path = str(self.paths["XML"][x].resolve())
        return self.xml_path

    def define_sbml_path(self):
        options = self.paths["SBML"]
        for i, item in enumerate(options):
            print(i, "||", item)
        x = int(input("Which model do you want? Start with 0: "))
        self.sbml_path = str(self.paths["SBML"][x].resolve())
        return self.sbml_path

    def load_model(self, path):
        self.model = read_sbml_model(path)
        return self.model

    def save_model(self):
        self.model = write_sbml_model()
