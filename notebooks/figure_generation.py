import os
import re
from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import cv2 as cv

pd.set_option("display.precision", 3)
pd.set_option("display.max_columns", 400)
pd.set_option("display.max_rows", 400)

import set_cwd
from src import Reader, Cleaner, BEAT, BNT, BAA

########################


def get_date_folder(target_date: str):

    tgt_path = os.path.join(os.getcwd(), "reports", "data", target_date)

    return tgt_path


def generate_pairgrid_dict(target_date: str):

    date_fdr = get_date_folder(target_date=target_date)

    pg_path_dict = {}
    strain_list = []

    for name in os.listdir(date_fdr):

        x = re.search(r"23-[0-9][0-9]_[0-9]+.[0-9]", name)

        strain_list.append(x.group())

    strain_set = set(strain_list)

    for strain in strain_set:

        strain_paths = []

        for pg in os.listdir(date_fdr):

            if strain in pg:
                strain_paths.append(pg)

        pg_path_dict[strain] = strain_paths

    return pg_path_dict


def generate_pair_grid_collage(target_date: str):

    strains_pg_dict = generate_pairgrid_dict(target_date=target_date)
    date_fdr = get_date_folder(target_date=target_date)

    collages_path = os.path.join(date_fdr, "collages")

    if not os.path.exists(collages_path):
        os.makedirs(collages_path)

    for st_tp_pgs in strains_pg_dict.keys():

        pth_list = []
        fig = plt.figure(
            figsize=(20, 20),
            layout="constrained",
        )

        for pg in strains_pg_dict[st_tp_pgs]:

            pg_path = os.path.join(date_fdr, pg)

            pth_list.append(pg_path)

        img1 = cv.imread(pth_list[0])
        img2 = cv.imread(pth_list[1])
        img3 = cv.imread(pth_list[2])
        img4 = cv.imread(pth_list[3])

        img1 = cv.cvtColor(img1, cv.COLOR_BGR2RGB)
        img2 = cv.cvtColor(img2, cv.COLOR_BGR2RGB)
        img3 = cv.cvtColor(img3, cv.COLOR_BGR2RGB)
        img4 = cv.cvtColor(img4, cv.COLOR_BGR2RGB)

        plt.subplot(2, 2, 1)
        plt.imshow(img1)
        plt.axis("off")

        plt.subplot(2, 2, 2)
        plt.imshow(img2)
        plt.axis("off")

        plt.subplot(2, 2, 3)
        plt.imshow(img3)
        plt.axis("off")

        plt.subplot(2, 2, 4)
        plt.imshow(img4)
        plt.axis("off")

        # * Save Image
        collage_path = os.path.join(collages_path, f"{st_tp_pgs}.png")
        fig.savefig(collage_path)

        plt.show()


x = generate_pair_grid_collage(target_date="2025_01_20")


########################
