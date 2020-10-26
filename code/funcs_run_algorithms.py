""" 
TO RUN PHYLONET
Create phylonet file with specific algorithm and settings
Create nexus files with convert_format.py
java -jar Phylonet_3.8.2.jar ./directory/file

TO RUN CASS
Before first time compile with:
javac ICass.java 
Then use
java ICass.java ./directory/file.tree

TO RUN TREE-CHILD NETWORKS
Use cargo to compile, then use
cargo run ./directory/inputfile -n -o ./directory/outputfile

TO RUN TRILONET
Change directory and filenames to all lowercase
java -jar TriLoNet.jar inputfile outputfile
"""

import subprocess
from funcs_format_conversions import *

datapath = "../data/"
resultspath = "../results/"
phylonetpath = "../ubuntu_software/phylonet/"
casspath = "../ubuntu_software/cass/"
dendropath = "../ubuntu_software/dendroscope/"
treechildpath = "../ubuntu_software/tree_child_code-master/"
temporalpath = "../ubuntu_software/temporal_hybridization_number-public_repo/"
trilonetpath = "../ubuntu_software/TriLoNet/TriLoNet/"
tril2netpath = "../ubuntu_software/TriL2Net-master/"

def run_phylonet(treefilename, phylonetfilename):
	runfilename = phylonetfilename.replace(".nex", "_") + treefilename.replace(".nex", ".txt")
	
	cmd1 = f"cd {phylonetpath}"
	cmd2 = f"cat ../{datapath}{treefilename} ./{phylonetfilename} > ./{runfilename}"
	cmd3 = f"java -jar PhyloNet_3.8.2.jar {runfilename}"
	
	print("Running PhyloNet")
	
	try:
		logfilename = treefilename.replace(".nex",".log")
		p1 = subprocess.run(f"{cmd1} && {cmd2} && {cmd3}", shell=True, capture_output=True)
		with open(f"../results/phylonet/{logfilename}", "w") as logfile:
			logfile.write(p1.stdout.decode("utf8"))
		phylonet_mpl_to_newick(logfilename)
		
	except:
		print("PhyloNet failed")

	finally:	
		pass
	
def run_cass(treefilename, maxtime):
	cmd1 = f"cd {casspath}"
	cmd2 = f"java ICass ../{datapath}{treefilename}"
	
	print("Running ICass")
	
	try:
		p1 = subprocess.run(f"{cmd1} && {cmd2}", shell=True, 
		capture_output=True, timeout=maxtime)
		print(p1.stderr.decode("utf8"))
		log = p1.stdout.decode("utf8")
		subprocess.run(f"mv {datapath}{treefilename.replace('.tree', '-network.tree')} ../results/cass/", shell=True)

	except subprocess.TimeoutExpired as timeout:
		log = timeout.stdout.decode("utf8")
		print("Timeout")
		
	finally:
		with open(f"../results/cass/{treefilename.replace('.dated.tree','.log')}", "w") as logfile:
			logfile.write(log)
		pass
		
def run_treechild(treefilename, maxtime):
	cmd1 = f"cd {treechildpath}code/tree_child/"
	cmd2 = f"cargo run ../../../{datapath}{treefilename} -n -o ../../../../results/treechild/{treefilename.replace('.tree','.net')}"
	
	print("Running Tree-Child Networks")
	
	try:
		p1 = subprocess.run(f"{cmd1} && {cmd2}", shell=True, timeout=maxtime)
		
	except subprocess.TimeoutExpired as timeout:
		print("Timeout")
		
	finally:
		pass

def run_temporal(treefilename, maxtime):
	cmd1 = f"cd {temporalpath}"
	cmd2 = f"./cherrypick_cpp -v ../{datapath}{treefilename}" # temporal
	cmd3 = f"./cherrypick_cpp -m 1 -v ../{datapath}{treefilename}" # semi-temporal
	
	print("Running Temporal Hybridization Number")

	p1 = subprocess.run(f"{cmd1} && {cmd2}", shell=True, capture_output=True)
	with open(f"../results/temporal/{treefilename.replace('.tree', '.log')}", "w") as logfile:
		logfile.write(p1.stdout.decode("utf8") + "\n")
	
	try:
		p2 = subprocess.run(f"{cmd1} && {cmd3}", shell=True, capture_output=True, timeout=maxtime)
		log = p2.stdout.decode("utf8")
		
	except subprocess.TimeoutExpired as timeout:
		log = timeout.stdout.decode("utf8")
		print("Timeout")

	finally:
		with open(f"../results/temporal/{treefilename.replace('.tree','.log')}", "a") as logfile:
			logfile.write(log)
		pass

def run_trilonet(infilename):
	outfilename = infilename.replace(".nex",".txt")
	
	cmd1 = f"cd {trilonetpath}"
	cmd2 = f"java -jar TriLoNet.jar /home/rosanne/Documents/minor_research_project/data/{infilename} /home/rosanne/Documents/minor_research_project/results/trilonet/{outfilename}"
	
	p1 = subprocess.run(f"{cmd1} && {cmd2}", shell=True)
	
	return(outfilename)
	
def run_tril2net(infilename):
	""" Run TriL2Net using trinets extracted from TriLoNet
	output as input """
	
	outfilename = infilename.replace(".tnet", "")
	
	cmd1 = f"cd {tril2netpath}"
	cmd2 = f"python3 tnet.py ../{resultspath}trilonet/{infilename} ../{resultspath}tril2net/{outfilename}"
	
	p1 = subprocess.run(f"{cmd1} && {cmd2}", shell=True)
	
	return(outfilename)
	
