import urllib.request
import subprocess
import re

# Set path names
datapath = "../data/"
resultspath = "../results/"
phylonetpath = "../ubuntu_software/phylonet/output/"
casspath = "../ubuntu_software/cass/"
dendropath = "../ubuntu_software/dendroscope/"
treechildpath = "../ubuntu_software/tree_child_code-master/"

def newick_to_nexus(filename):
	""" Converts file with multiple trees in Newick format 
	(one tree per line) to NEXUS format (with numbered
	trees) as suitable for PhyloNet. """
	
	infile = open(datapath+filename,"r")
	outfile = open(datapath+filename.replace(".tree",".nex"),"w")
	outfile.write("#NEXUS\nBEGIN trees;\n")

	i=1
	for line in infile:
		outfile.write(f"TREE Tree{i} = {line}")
		i += 1
		
	outfile.write("\nEND;\n")
	infile.close()
	outfile.close()
	
def phylip_to_nexus(infilename):
	""" Converts multiple sequence alignment file in Phylip
	format to NEXUS format as suitable for TriLoNet """
	
	infile = open(datapath+infilename, "r")
	outfilename = infilename.replace(".ph", ".nex").lower() # filename should be lowercase for TriLoNet on Linux OS
	outfile = open(datapath+outfilename, "w")
	header = infile.readline().strip().split(" ")
	ntax = header[0]
	nchar = header[1]
	
	outfile.write(f"#NEXUS\n\nBEGIN taxa;\n\tDIMENSIONS NTAX={ntax};\n\tTAXLABELS\n")
	for line in infile:
		outfile.write(f"\t\t{line.split()[0]}\n")
		
	infile.close()
	infile = open(datapath+infilename, "r")
	infile.readline()
	
	outfile.write(f"\t;\nEND;\n\nBEGIN characters;\n\tDIMENSIONS NCHAR={nchar};\n\tFORMAT datatype=nucleotide gap=- mising=? matchchar=.;\n\tMATRIX\n")
	
	for line in infile:
		outfile.write(line)
		
	outfile.write("\t;\nEND;")
		
	infile.close()
	outfile.close()
	
	return(outfilename)
	
def phylonet_mpl_to_newick(filename):
	""" Converts PhyloNet Infer Network MPL output to
	multiple networks in Newick format as
	readable for Dendroscope. """
	
	infile = open(f"{resultspath}phylonet/{filename}","r")
	outfile = open(f"{resultspath}phylonet/{filename.replace('.log', '.net')}", "w")
	
	i=0
	for line in infile:
		if "Visualize in Dendroscope" in line:
			i += 1
			network = line.split(" : ")[1]
			outfile.write(network)

	infile.close()
	outfile.close()

def find_dates(filename, datefilename):
	""" Finds collection dates of taxon names in NEXUS
	multiple sequence alignment file (for improved
	CASS). Assumes that taxon names start with NCBI
	accession number. Retrieves collection data 
	from NCBI. """
	
	infile = open(datapath+filename,"r")
	outfile = open(datapath+datefilename,"w")
	
	taxlabels = False
	while taxlabels == False:
		line = infile.readline()
		if "ntax" in line.lower():
			ntax = int(line.split("=")[1].strip()[:-1])
		if "taxlabels" in line.lower():
			taxlabels = True
	
	for t in range(0, ntax):
		taxon = infile.readline().strip().strip(";")
		if taxon.startswith("NC"):
			acc = taxon
		else:
			acc = taxon.split("_")[0]

		url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=native&retmode=xml"
		fh = urllib.request.urlopen(url)
		info = fh.read().decode("utf8").split("\n")
		fh.close()
		dateFound = False
		for i in range(0, len(info)):
			line = info[i]
			if "collection-date" in line:
				dateFound = True
				break
		if dateFound == True:
			date = info[i+1].split(">")[1].split("<")[0]
		else:
			date = "0000"
			print(f"collection date for accession number {acc} cannot be found")
		outfile.write(f"{taxon} {date} \n")
	
	infile.close()
	outfile.close()

def add_dates(infilename, datefilename):
	infile = open(datapath+infilename, "r")
	with open(datapath+datefilename, "r") as datefile:
		taxadates = datefile.read().strip().split("\n")
	outfile = open(datapath+infilename.replace(".tree", ".dated.tree"), "w")
	
	for line in infile:
		for taxon in taxadates:
			name = taxon.split()[0]
			date = taxon.split()[1]
			if name in line:
				line = line.replace(name, name+"/"+date)

		outfile.write(line)
		
	infile.close()
	datefile.close()
	outfile.close()
	
def replace_scientific_numbers(infilename):
	""" Replaces scientific numbers in Newick trees (converted
	by Dendroscope for very small branch lengths) with 0.0 """
	# TO DO maybe convert to actual non-scientific number
	# but not necessary for the algorithms I use
	
	infile = open(datapath+infilename, "r")
	outfile = open(datapath + "a_" + infilename, "w")
	
	for line in infile:
		newline = line
		line = re.split(":|,|\)", line) # split on ":" and "," and ")"
		for word in line:
			if "E-" in word:
				newline = newline.replace(word, "0.0")		
		outfile.write(newline)
		
	infile.close()
	outfile.close()
	
def remove_internal_node_labels(infilename):
	""" Removes internal node labels such as bootstrap support values
	while retaining edge labels such as branch length
	from a file in newick format (one tree per line) """
	
	infile = open(datapath+infilename, "r")
	outfile = open(datapath + infilename.replace(".tree", "_noBS.tree"), "w")
	
	for line in infile:
		line = line.strip().split(",")
		newline = ""
		for word in line:
			word = word.split(":")
			newword = ""
			for letter in word:
				close_bracket = letter.count(")")
				newletter = letter.split(")")[0]
				newword = newword + newletter + close_bracket*")" + ":"
			newline = newline + newword.strip(":") + ","
		outfile.write(newline.strip(",") + ";\n")
		
	infile.close()
	outfile.close()

def remove_edge_labels(infilename):
	""" Removes all edge and internal node labels such as branch 
	length (edge label) and bootstrap support values
	from a file in newick format (one tree per line) """
	
	infile = open(datapath+infilename, "r")
	outfile = open(datapath + infilename.replace(".tree", "_nolen.tree"), "w")
	
	for line in infile:
		line = line.strip().split(",")
		newline = ""
		for word in line:
			leaf_label = word.split(":")[0]
			close_bracket = word.count(")")
			newline = newline + leaf_label + close_bracket*")" + ","
		outfile.write(newline.strip(",") + ";\n")
		
	infile.close()
	outfile.close()

def midpoint_root(infilename):
	infile = open(datapath+infilename, "r")
	outfilename = infilename.replace(".tree", "_rooted.tree")
	
	with open(dendropath+"temp_cmdfile.txt", "w") as cmdfile:
		cmdfile.write(f"open file='/home/rosanne/Documents/minor_research_project/{datapath.replace('../','')}{infilename}';\n")
		cmdfile.write("midpointroot;\n")
		cmdfile.write(f"save format=newick file='/home/rosanne/Documents/minor_research_project/{datapath.replace('../','')}{outfilename}';\nquit")
		
	cmd1 = "cd " + dendropath
	cmd2 = "./Dendroscope -g -c ./temp_cmdfile.txt"
	
	p1 = subprocess.call(cmd1 + " && " + cmd2, shell=True)

def contract_edges(infilename, len_threshold=None, supp_threshold=None):
	infile = open(datapath+infilename, "r")
	outfilename = infilename.replace("_branchlen_BSsupport", "")

	with open(dendropath+"temp_cmdfile.txt", "w") as cmdfile:
		cmdfile.write(f"open file='/home/rosanne/Documents/minor_research_project/{datapath.replace('../','')}{infilename}';\n")
		if len_threshold is not None:
			cmdfile.write(f"select edges=short threshold={len_threshold};\n")
			cmdfile.write(f"remove edges=selected;\n")
			outfilename = outfilename.replace(".tree", f"_len{len_threshold}.tree")

		if supp_threshold is not None:
			cmdfile.write(f"set interpretInternalNodeLabelsAsEdgeLabels=true;\n")
			cmdfile.write(f"contractEdges minSupport={supp_threshold};\n")
			outfilename = outfilename.replace(".tree", f"_supp{supp_threshold}.tree")
		
		cmdfile.write(f"save format=newick file='/home/rosanne/Documents/minor_research_project/{datapath.replace('../','')}{outfilename}';\nquit")
		
	cmd1 = "cd " + dendropath
	cmd2 = "./Dendroscope -g -c ./temp_cmdfile.txt"
	
	p1 = subprocess.call(cmd1 + " && " + cmd2, shell=True)
	
	return(outfilename)
	
def nonbinary_to_binary(infilename):
	with open(datapath+infilename, "r") as infile:
		line = infile.readline().split(",")
		nleafs = len(line)
	
	cmd1 = f"cd {treechildpath}code/data_gen/real_world"
	cmd2 = f"MakeTestData-exe {nleafs} ../../../../{datapath}/{infilename} ../../../../{datapath}/{infilename.replace('.tree', '_binary.tree')}"

	p1 = subprocess.run(cmd1 + " && " + cmd2, shell=True, capture_output=True)
	
	logfilename = infilename.replace(".tree", ".log")
	with open(f"{resultspath}nonbinary-to-binary/{logfilename}", "w") as logfile:
		logfile.write(p1.stdout.decode("utf8"))

	# remove quotation marks from output file
	with open(f"{datapath}{infilename.replace('.tree', '_binary.tree')}", "r") as outfile:
		out = outfile.read().replace('"','')

	with open(f"{datapath}{infilename.replace('.tree', '_binary.tree')}", "w") as outfile:
		outfile.write(out)
		
def trilonet_output_to_network(infilename):
	""" Creates a .net file with the eNewick network string
	as found in the .txt output file generated by TriLoNet """
	
	infile = open(f"{resultspath}trilonet/{infilename}", "r")
	
	output = False
	for line in infile:
		if output:
			network = line.strip()
			break
		if "Output eNewick String Short" in line:
			output = True
			
	infile.close()
	with open(f"{resultspath}trilonet/{infilename.replace('.txt', '.net')}", "w") as outfile:
		outfile.write(network)
		
def trilonet_output_to_trinets(infilename):
	""" Creates a .tnet file with trinets from the .txt
	output file generated by TriLoNet """

	infile = open(f"{resultspath}trilonet/{infilename}", "r")
	outfilename = infilename.replace(".txt", ".tnet")
	outfile = open(f"{resultspath}trilonet/{outfilename}", "w")
	
	for line in infile:
		if line.startswith("Tr"):
			outfile.write(line)
			
	infile.close()
	outfile.close()
	
	return(outfilename)
	
	

