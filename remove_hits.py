import pandas as pd
import os
import argparse

def remove_not_satellite(path):
    """
    :param path: path to RE2 archive folder
    :return: single FILTERED_SATELLITES.fasta file containing all found satellites
    """
    tarean_report_path = os.path.join(path, "tarean_report.html")
    fasta_1_path = os.path.join(path, "TAREAN_consensus_rank_1.fasta")
    fasta_2_path = os.path.join(path, "TAREAN_consensus_rank_2.fasta")
    is_fasta_1 = None
    is_fasta_2 = None
    with open(fasta_1_path) as file:
        fasta_content = file.read().strip()
        if fasta_content:
            is_fasta_1 = 1

    with open(fasta_2_path) as file:
        fasta_content = file.read().strip()
        if fasta_content:
            if is_fasta_1:
                is_fasta_2 = 3
            else:
                is_fasta_2 = 1

    hc_lc = [is_fasta_1, is_fasta_2]
    if os.path.exists(os.path.join(path, "FILTERED_SATELLITES.fasta")):
        print("File FILTERED_SATELLITES.fasta already exists! Remove before proceeding")
        return 0
    tarean_report = pd.read_html(tarean_report_path)
    for flag in hc_lc:
        try:
            table = tarean_report[flag]
        except TypeError:
            continue
        for row in table.iterrows():
            hit = row[1]["Similarity hits [above 0.1%]"]
            if isinstance(hit, str):
                continue
            cluster_n = row[1]["Cluster"]
            cluster_len = row[1]["Consensus length"]
            cluster_seq = row[1]["Consensus"]
            with open(os.path.join(path, "FILTERED_SATELLITES.fasta"), "a") as outfile:
                outfile.write(">CL" + str(cluster_n) + "_" + str(cluster_len) + "\n")
                outfile.write(cluster_seq)
                outfile.write("\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="Path to RE2 run folder")

    args = parser.parse_args()

    remove_not_satellite(args.path)

if __name__ == "__main__":
    main()
