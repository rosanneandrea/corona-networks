import subprocess
import urllib.request
import argparse
import re
from funcs_format_conversions import *
from funcs_run_algorithms import *

def main(tree_fname, tree_bs_fname, date_fname, phylonet_fname, len_threshold, supp_threshold, treechild, temporal, icass):
	maxtime = 5*60
	switch = False
		
	if len_threshold or supp_threshold:
		tree_bs_fname = contract_edges(tree_bs_fname, len_threshold, supp_threshold)
		switch = True

	if icass:
		remove_internal_node_labels(tree_bs_fname)
		add_dates(tree_bs_fname.replace(".tree", "_noBS.tree"), date_fname)
		run_cass(tree_bs_fname.replace(".tree", "_noBS.dated.tree"), maxtime)
		subprocess.run("killall -9 java", shell=True)
		
	if treechild or temporal:
		if switch:
			nonbinary_to_binary(tree_bs_fname)
			tree_fname = tree_bs_fname.replace(".tree", "_binary.tree")
		else:
			remove_edge_labels(tree_bs_fname)
			tree_fname = tree_bs_fname.replace(".tree", "_nolen.tree")
	
	if treechild:
		run_treechild(tree_fname, maxtime)
		subprocess.run("killall -9 tc_seq", shell=True)
		
	if temporal:
		run_temporal(tree_fname, maxtime)		
		
	if phylonet_fname:
		phylonet_fname = "MPL.nex"
		newick_to_nexus(tree_bs_fname)
		run_phylonet(tree_bs_fname.replace(".tree", ".nex"), phylonet_fname)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Preprocess tree files and construct phylogenetic networks with various algorithms.")
	parser.add_argument("-trees", "--trees_plain", help="file with Newick trees including branch lengths", required=False)
	parser.add_argument("-bs", "--trees_bootstrap", help="file with Newick trees including branch lengths and bootstrap support values", required=True)
	parser.add_argument("-dates", "--collection_dates", help="file containing on each line sequence ID and its collection year (space-separated), required for iCASS", required=False)
	parser.add_argument("-len", "--branch_length", help="threshold for edge contraction based on branch length", required=False)
	parser.add_argument("-supp", "--bootstrap_support", help="threshold for edge contraction based on branch length", required=False)
	parser.add_argument("-tcn", "--treechild", dest="treechild", action="store_true", help="flag to run Tree-Child Networks algorithm")
	parser.add_argument("-temp", "--temporal", dest="temporal", action="store_true", help="flag to run Temporal Hybridization Number algorithm")
	parser.add_argument("-icass", "--icass", dest="icass", action="store_true", help="flag to run improved CASS algorithm")
	parser.add_argument("-phylo", "--phylonet", help="PhyloNet filename, also serves as a flag for this algorithm")
	
	args = parser.parse_args()
	
	# Input files
	tree_fname = args.trees_plain
	tree_bs_fname = args.trees_bootstrap
	date_fname = args.collection_dates
	phylonet_fname = args.phylonet # also serves as a flag for PhyloNet
	
	# Thresholds for edge contraction
	len_threshold = args.branch_length
	supp_threshold = args.bootstrap_support

	# Algorithm flags
	treechild = args.treechild
	temporal = args.temporal
	icass = args.icass

	main(tree_fname, tree_bs_fname, date_fname, phylonet_fname, len_threshold, supp_threshold, treechild, temporal, icass)

