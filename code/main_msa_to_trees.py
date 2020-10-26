import subprocess
import urllib.request
import argparse
import re
from funcs_format_conversions import *
from funcs_msa_to_trees import *

# Set path names
datapath = "../data/"
phylonetpath = "../ubuntu_software/phylonet/"
casspath = "../ubuntu_software/cass/"
dendropath = "../ubuntu_software/dendroscope/"
treechildpath = "../ubuntu_software/tree_child_code-master/"
resultspath = "/home/rosanne/Documents/minor_research_project/results/"
raxmlpath = "../ubuntu_software/standard-RAxML-master/"

"""
1) Remove sequences from Phylip MSA
2) Create RAxML exclude files from gene locations or breakpoints
3) Build trees with RAxML
4) Midpoint rooting both output best trees (with and without BS support values)
"""

def main(version, msaphylip_fname, msanexus_fname, excl_seqs_fname, bp_fname):

	if excl_seqs_fname:
		exclude_sequences(msaphylip_fname, excl_seqs_fname)
		msaphylip_fname = msaphylip_fname.split(".")[0] + "_" + excl_seqs_fname.replace(".txt", ".ph")
		version = version + "_" + excl_seqs_fname.split(".")[0]
	
	breakpoints = get_breakpoints(bp_fname)
	
	create_excludefiles(version, msaphylip_fname, breakpoints)
	run_raxml(version, msaphylip_fname, breakpoints)
	cat_raxml_trees(version, bp_fname)

	tree_fname = version + ".tree"
	tree_bs_fname = version + "_BS.tree"
	midpoint_root(tree_fname)
	midpoint_root(tree_bs_fname)
	replace_scientific_numbers(tree_fname.replace(".tree","_rooted.tree"))
	replace_scientific_numbers(tree_bs_fname.replace(".tree","_rooted.tree"))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Build maximum-likelihood trees for different genetic locations from multiple sequence alignment file in Phylip format. With option of excluding sequences.")
	parser.add_argument("-v", "--version", help="version name", required=True)
	parser.add_argument("-msap", "--msa_phylip", help="multiple sequence alignment file in Phylip format", required=True)
	parser.add_argument("-excl", "--exclude_seqs", help="file with IDs of sequences to exclude before building trees", required=False)
	parser.add_argument("-bp", "--breakpoint_locations", help="file with breakpoint locations (e.g. genes) to build separate trees for, on each line: name begin end (space-separated)", required=False)
	
	args = parser.parse_args()
	
	# Input files
	msaphylip_fname = args.msa_phylip
	excl_seqs_fname = args.exclude_seqs
	bp_fname = args.breakpoint_locations
	
	version = args.version

	main(version, msaphylip_fname, msanexus_fname, excl_seqs_fname, bp_fname)

