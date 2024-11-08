from Bio import SeqIO
import os
import argparse

def make_fasta_file(path: str, to_save: str, n: str):
    """
    :param path: path to RE2 run_folder
    :param to_save: where to save resulted file
    :param n: n-times to concatenate satellite sequences
    :return: fasta file with
    """
    tarean_files = [file for file in os.listdir(path) if "TAREAN_consensus" in file]
    n_mer = n + "_mer.fasta"
    with open(f"{os.path.join(to_save, n_mer)}", "w") as outfile:
        for tarean_file in tarean_files:
            path_to_tarean_file = os.path.join(path, tarean_file)
            parser = SeqIO.parse(path_to_tarean_file, "fasta")

            for entry in parser:
                outfile.write(">" + str(entry.id) + "\n")
                outfile.write(str(entry.seq)*int(n))
                outfile.write("\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="path to RE2 run folder")
    parser.add_argument("--to_save", required=True, help="path to resulted fasta file")
    parser.add_argument("--n_mers", default="5", help="how many satellites to concatenate. Default 5")
    args = parser.parse_args()

    make_fasta_file(args.path, args.to_save, args.n_mers)

if __name__ == "__main__":
    main()
