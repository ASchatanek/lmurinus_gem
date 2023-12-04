from pathlib import Path
import pandas as pd


class Reader:
    def __init__(self) -> None:
        self.find_folder()

    def find_folder(self) -> Path:
        self.biolog_folder = Path.cwd() / "data" / "biolog_results"
        self.biolog_folder.resolve()

        return self.biolog_folder

    def find_file(self, file_name: str) -> Path:
        self.file_path = self.biolog_folder / file_name
        self.file_path.resolve()
        return self.file_path

    def find_dataframe(self, file_name: str) -> pd.DataFrame:
        self.find_file(file_name=file_name)

        file_sheet_names = pd.ExcelFile(self.file_path).sheet_names

        print("Available Sheet Names:")
        for index, sheet_name in enumerate(file_sheet_names):
            print(f"{index}. {sheet_name}")
            print("")

        target_sheet = input("Which sheet should be loaded?")

        self.df = pd.read_excel(self.file_path, target_sheet)

        found_df = self.df.copy()

        return found_df

    def get_compounds_list(self, compounds_df) -> list:
        self.compounds = list(compounds_df["Compound"])
        return self.compounds
