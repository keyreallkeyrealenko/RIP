from Bio import SeqIO
import subprocess
import argparse
import os
from create_dotplots import find_clusters

def assemble_reads(path, n_threads):
    """
    :param path: path to a RE2 run folder
    :param n_threads: number of threads for assembly
    :return: assembled clusters possibly containing satellites
    """
    clusters_to_assemble = find_clusters(path)
    path_to_clusters = os.path.join(path,"seqclust", "clustering", "clusters")
    for cluster_dir in clusters_to_assemble:
        current_cl_path = os.path.join(path_to_clusters, cluster_dir)
        if os.path.exists(os.path.join(current_cl_path, "assembly")):
            print(f"{current_cl_path} Already exist, remove before proceeding")
            continue
        os.makedirs(os.path.join(current_cl_path, "assembly"))

        input_file = os.path.join(current_cl_path, "reads.fasta")
        forward_file = os.path.join(current_cl_path, "assembly", "forward_reads.fasta")
        reverse_file = os.path.join(current_cl_path, "assembly", "reverse_reads.fasta")
        with open(forward_file, "w") as fwd, open(reverse_file, "w") as rev:
            parser = list(SeqIO.parse(input_file, "fasta"))
            for idx in range(len(parser) - 1):
                if parser[idx].id.replace("f", "").replace("r", "") == parser[idx + 1].id.replace("f", "").replace("r",
                                                                                                                   ""):
                    if parser[idx].id.endswith("f"):
                        SeqIO.write(parser[idx], fwd, "fasta")
                        SeqIO.write(parser[idx + 1], rev, "fasta")
                    else:
                        SeqIO.write(parser[idx + 1], fwd, "fasta")
                        SeqIO.write(parser[idx], rev, "fasta")

        megahit_command = ["megahit", "-1" , os.path.join(current_cl_path, "assembly", "forward_reads.fasta"), "-2",
                           os.path.join(current_cl_path, "assembly", "reverse_reads.fasta"), "-o",
                           os.path.join(current_cl_path, "assembly", "megahit"), "--presets", "meta-large", "-t",
                           n_threads]
        print(megahit_command)
        subprocess.run(megahit_command)

        spades_command = ["spades.py", "--isolate", "--threads", n_threads, "--cov-cutoff", "auto", "-1",
                          os.path.join(current_cl_path, "assembly", "forward_reads.fasta"), "-2",
                          os.path.join(current_cl_path, "assembly", "reverse_reads.fasta"), "--trusted-contigs",
                          os.path.join(current_cl_path, "assembly", "megahit", "final.contigs.fa"), "-o",
                          os.path.join(current_cl_path, "assembly")]
        subprocess.run(spades_command)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="Path to the run folder")
    parser.add_argument("--n_threads", default="2", help="Number of threads for assembly")

    args = parser.parse_args()
    assemble_reads(args.path, args.n_threads)

if __name__ == "__main__":
    main()
