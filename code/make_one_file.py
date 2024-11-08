import os
import argparse

def concatenate_files(path: str):
    run_folders = [folder for folder in os.listdir(path) if "run" in folder]
    if os.path.exists(os.path.join(path, "Satellites.fasta")):
        print(f"Satellites.fasta already exists in {path}. Remove before proceeding.")
        return 0
    for folder in run_folders:
        dir_name = os.path.basename(os.path.dirname(os.path.join(path, folder)))
        path_to_file = os.path.join(path, folder, "FILTERED_SATELLITES.fasta")
        if not os.path.exists(path_to_file):
            print(f"FILTERED_SATELLITES.fasta in {os.path.join(path, folder)} doesn't exist. Leaving.")
            return 0
        with open(os.path.join(path, "Satellites.fasta"), "a") as outfile:
            with open(os.path.join(path, folder, "FILTERED_SATELLITES.fasta")) as infile:
                for line in infile.readlines():
                    if line.startswith(">"):
                        line = ">" + dir_name + "_" + line.strip().replace(">", "") + "_" + folder + "\n"
                    outfile.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="path to the root folder")
    args = parser.parse_args()
    concatenate_files(args.path)

if __name__ == "__main__":
    main()
