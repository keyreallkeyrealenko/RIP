# <center> RIP </center>
### RepeatExplorer2 Inspector Program

Here, you can find a collection of scripts to parse your RepeatExplorer2 output, particularly for putative satellites.

As a result, RE2 provides a folder containing putative repetitive DNA, including satellites and transposons. However, it is not trivial to validate these findings. WWe present an approach based on several methods:

- Visualize satellite coverage (`make_n_mer_fasta.py` and `visualize_n_mers.py` scripts). Satellites have a homogenized nature,
so if you align reads to concatenated satellite sequences (e.g., 5 monomers concatenated into one sequence), the coverage across the entire n-mer should be even. 
- Assemble reads assigned to a cluster by RE2 and create a dot-plot using the most covered and the longest scaffolds 
(`assemble_clusters.py`, `make_one_file` and `create_dotplots` scripts). 