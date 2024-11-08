import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

def visualize_n_mers(bed_file: str, to_save: str):
    """
    :param bed_file: file
    :param to_save:
    :return:
    """
    df = pd.read_csv(bed_file, sep="\t", header=None, names=["Cl", "pos", "cov"])

    n_unique_clusters = df["Cl"].nunique()
    fig, axes = plt.subplots(ncols=1, nrows=n_unique_clusters, figsize=(6, 3*n_unique_clusters))
    unique_clusters = list(df["Cl"].unique())

    for ax in axes:
        cl = unique_clusters.pop()
        m_len = int(cl.split("_")[-1].replace("nt", ""))
        subset = df[df["Cl"] == cl]

        ax.fill_between(subset["pos"], subset["cov"], alpha=0.8, color="#99a3a4")
        ax.plot(subset["pos"], subset["cov"], color='black', linewidth=0.5)
        x_unique = np.unique(subset["pos"])
        num_xticks = len(x_unique) // 10
        ax.set_xticks(ticks=x_unique[::num_xticks])
        ax.set_title(subset["Cl"].unique()[0])
        ax.set_ylabel("")
        ax.set_xlabel("")

        v_lines_height = subset["cov"].max()
        for i in range(1, 5):
            ax.vlines(x=m_len*i, ymin=0, ymax=v_lines_height, color = 'r')


    plt.savefig(f"{os.path.join(to_save, "5_merPLOT.png")}", dpi=150)
    plt.subplots_adjust(hspace=0.5)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed_file", required=True, help="path to bed file")
    parser.add_argument("--to_save", required=True, help="path to save a plot")
    args = parser.parse_args()

    visualize_n_mers(args.bed_file, args.to_save)

if __name__ == "__main__":
    main()
