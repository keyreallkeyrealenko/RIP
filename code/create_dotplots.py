from Bio import SeqIO
import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
from Levenshtein import distance as lev_distance
import subprocess
from remove_hits import remove_not_satellite

def find_clusters(path: str):
    """
    path - path to RE2 run folder
    Finds only those clusters to work with. They're in the TAREAN_consensus* files
    """
    if os.path.exists(os.path.join(path, "FILTERED_SATELLITES.fasta")):
        sat_path = os.path.join(path, "FILTERED_SATELLITES.fasta")
    else:
        remove_not_satellite(path)
        sat_path = os.path.join(path, "FILTERED_SATELLITES.fasta")

    clusters = []

    with open(os.path.join(path, sat_path)) as sat_file:
        for line in sat_file.readlines():
            if line.startswith(">"):
                cl = line.split("_")[0].replace(">CL", "")
                cl = "dir_CL" + cl
                while len(cl) < 10:
                    tmp_list = cl.split("CL")
                    cl = tmp_list[0] + "CL" + "0" + tmp_list[1]

                clusters.append(cl)

    return clusters

def find_best_scaffolds(path: str, clusters: list):
    """
    :param path: path to RE2 run folder
    :param clusters: clusters to works with (fasta headers)
    :return: a dictionary with most covered scaffold and the longest scaffold
    """
    path_to_clusters = os.path.join(path, "seqclust", "clustering", "clusters")

    most_covered_scaffolds = {}
    for cluster_dir in clusters:
        current_cl_path = os.path.join(path_to_clusters, cluster_dir, "assembly")
        try:
            parser = SeqIO.parse(os.path.join(current_cl_path, "scaffolds.fasta"), "fasta")

            cov_threshold = 0
            len_threshold = 0
            len_seq = None
            cov_seq = None
            for entry in parser:
                cur_len_threshold = float(entry.id.split("_")[-3])
                cur_cov_threshold = float(entry.id.split("_")[-1])
                if cur_len_threshold >= len_threshold:
                    len_threshold = cur_len_threshold
                    len_seq = str(entry.seq)
                if cur_cov_threshold >= cov_threshold:
                    cov_threshold = cur_cov_threshold
                    cov_seq = str(entry.seq)
            most_covered_scaffolds[cluster_dir] = [len_seq, cov_seq]
        except FileNotFoundError:
            most_covered_scaffolds[cluster_dir] = None

    return most_covered_scaffolds


def create_dotplot(path: str, clusters: list, most_cov_scaffolds: dict,  window: int = 10, max_distance: int = 1):
    """
    :param path: path to RE2 run folder
    :param clusters: clusters to work with.
    :param most_cov_scaffolds: a dictionary with the most covered and the longest scaffolds
    :param window: slicing window to count distance (Levenstein algorithm). default 10
    :param max_distance: maximum distance between two slices. 1 means 1 alteration in seqs. Default 1
    :return: creates dotplots from the most covered and the longest scaffold for each assembly
    """
    path_to_clusters = os.path.join(path, "seqclust", "clustering", "clusters")
    for cluster_dir in clusters:

        if not most_cov_scaffolds[cluster_dir]:
            continue

        current_cl_path = os.path.join(path_to_clusters, cluster_dir, "assembly")
        if os.path.exists((os.path.join(current_cl_path, "DOTPLOT"))):
            print(f"a DOTPLOT folder already exists in a {cluster_dir}! Remove before proceeding")
            return 0
        os.makedirs(os.path.join(current_cl_path, "DOTPLOT"))
        cons_dir = os.path.join(current_cl_path, "DOTPLOT")
        len_seq = most_cov_scaffolds[cluster_dir][0]
        cov_seq = most_cov_scaffolds[cluster_dir][1]

        # length seq:
        length = len(len_seq)
        matrix = np.zeros((length, length))

        for i in range(length - window + 1):
            for j in range(length - window + 1):
                if lev_distance(len_seq[i:i + window], len_seq[j:j + window]) <= max_distance:
                    matrix[i:i + window, j:j + window] = 1

        plt.imshow(matrix, cmap='Greys', interpolation='none')
        plt.yticks(range(1, length, length // 25), size=8)
        plt.xticks(range(1, length, length // 25), rotation=90, size=8)
        plt.savefig(os.path.join(cons_dir, "dotplot_longest_scaffold.png"), dpi=150)

        # coverage seq:
        length = len(cov_seq)
        matrix = np.zeros((length, length))

        for i in range(length - window + 1):
            for j in range(length - window + 1):
                if lev_distance(cov_seq[i:i + window], cov_seq[j:j + window]) <= max_distance:
                    matrix[i:i + window, j:j + window] = 1

        plt.imshow(matrix, cmap='Greys', interpolation='none')
        plt.yticks(range(1, length, length // 25), size=8)
        plt.xticks(range(1, length, length // 25), rotation=90, size=8)
        plt.savefig(os.path.join(cons_dir, "dotplot_most_covered_scaffold.png"), dpi=150)


def make_consensus(clusters: list, path: str, most_cov_scaffolds: dict):
    """
    clusters - clusters to work with.
    Creates consensus from most covered scaffold for each cluster. Uses TRF to make a consensus.
    """
    path_to_clusters = os.path.join(path, "seqclust", "clustering", "clusters")
    for cluster_dir in clusters:

        if not most_cov_scaffolds[cluster_dir]:
            continue

        current_cl_path = os.path.join(path_to_clusters, cluster_dir, "assembly")
        os.makedirs(os.path.join(current_cl_path, "CONSENSUS"))
        cons_dir = os.path.join(current_cl_path, "CONSENSUS")
        len_seq = most_cov_scaffolds[cluster_dir][0]
        cov_seq = most_cov_scaffolds[cluster_dir][1]

        with open(os.path.join(cons_dir, "longest_scaffold.fasta"), "w") as infile:
            infile.write(">NODE_1"+"\n")
            infile.write(len_seq)

        with open(os.path.join(cons_dir, "most_covered_scaffold.fasta"), "w") as infile:
            infile.write(">NODE_1"+"\n")
            infile.write(cov_seq)

        trf_command_1 = ["trf", os.path.join(cons_dir, "most_covered_scaffold.fasta"), "2", "7", "7", "80", "10", "1",
                         "5000", "-d"]
        trf_command_2 = ["trf", os.path.join(cons_dir, "longest_scaffold.fasta"), "2", "7", "7", "80", "10", "1",
                         "5000", "-d"]
        os.chdir(os.path.join(cons_dir))
        subprocess.run(trf_command_1)
        subprocess.run(trf_command_2)
        try:
            with open(os.path.join(cons_dir, "longest_scaffold.fasta.2.7.7.80.10.1.5000.dat")) as infile:
                seq_threshold = 0
                for line in infile.readlines():
                    line = line.split(" ")
                    if line[0].isdigit():
                        quality = int(line[7])
                        if quality > seq_threshold:
                            seq_threshold = quality
                            seq = line[-2]

            with open(os.path.join(cons_dir, "longest_CONSENSUS.fasta"), "w") as outfile:
                outfile.write(">Consensus_1"+"\n")
                outfile.write(seq)
        except:
            pass

        try:
            with open(os.path.join(cons_dir, "most_covered_scaffold.fasta.2.7.7.80.10.1.5000.dat")) as infile:
                seq_threshold = 0
                for line in infile.readlines():
                    line = line.split(" ")
                    if line[0].isdigit():
                        quality = int(line[7])
                        if quality > seq_threshold:
                            seq_threshold = quality
                            seq = line[-2]

            with open(os.path.join(cons_dir, "most_covered_CONSENSUS.fasta"), "w") as outfile:
                outfile.write(">Consensus_2" + "\n")
                outfile.write(seq)
        except:
            pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="Path to the run folder")
    parser.add_argument("--TRF", help="Create consensus using TRF. TRF must be in PATH variable",
                        action="store_true")
    args = parser.parse_args()

    clusters = find_clusters(args.path)
    best_scaffolds = find_best_scaffolds(path=args.path, clusters=clusters)
    create_dotplot(path=args.path, clusters=clusters, most_cov_scaffolds=best_scaffolds)
    if args.TRF:
        make_consensus(path=args.path, clusters=clusters, most_cov_scaffolds=best_scaffolds)

if __name__ == "__main__":
    main()
