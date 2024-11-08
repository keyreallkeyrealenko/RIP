# <center> RIP </center>
### RepeatExplorer2 Inspector Program

Here, you can find a collection of scripts to parse your RepeatExplorer2 output, particularly for putative satellites.

As a result, RE2 provides a folder containing putative repetitive DNA, including satellites and transposons. However, it is not trivial to validate these findings. We present an approach based on several approaches:

- Visualize satellite coverage (`make_n_mer_fasta.py` and `visualize_n_mers.py` scripts). Satellites have homogenized nature, 
so, if you align reads on concatenated satellites sequences (e.g. 5 monomers concatenated in one sequence), the coverage across 
the whole n-mer should be even. 
- 