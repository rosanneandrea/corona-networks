import subprocess
import argparse

datapath = "../data/"
resultspath = "../results/"
trilonetpath = "../ubuntu_software/TriLoNet/TriLoNet"
tril2netpath = "../ubuntu_software/TriL2Net-master/"

from funcs_format_conversions import trilonet_output_to_network, trilonet_output_to_trinets, phylip_to_nexus
from funcs_run_algorithms import run_trilonet, run_tril2net

def main(msa_fname):
	#if msap_fname:
		# convert msa from phylip to nexus
	#	msa_fname = phylip_to_nexus(msap_fname)
	msa_fname = msap_fname.replace(".ph", ".nex").lower()
	
	#out_fname = run_trilonet(msa_fname)
	#trilonet_output_to_network(out_fname)
	#trinet_fname = trilonet_output_to_trinets(out_fname)
	trinet_fname = msa_fname.replace(".nex", ".tnet")
	run_tril2net(trinet_fname)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Construct phylogenetic networks with TriLoNet and TriL2Net from multiple sequence alignment file in Nexus or Phylip format.")
	parser.add_argument("-msan", "--msa_nexus", help="multiple sequence alignment file in single-line Nexus format", required=False)
	parser.add_argument("-msap", "--msa_phylip", help="multiple sequence alignment file in Phylip format", required=False)
	
	args = parser.parse_args()
	msa_fname = args.msa_nexus
	msap_fname = args.msa_phylip
	
	main(msa_fname)
