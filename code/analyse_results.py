import subprocess
import re

resultspath = "../results/"
datapath = "../data/"
dendropath = "../ubuntu_software/dendroscope/"

def count_reticulation_number(infilename):
	# file containing networks in extended
	# newick format, one per line
	infile = open(infilename, "r")
	retnumlist = []
	
	for line in infile:
		line = line.strip().split("#")[1:] # skip part before first #
		reticulations = {}
		
		for word in line:
			retic = word[1]
			try:
				# add one count to this reticulation
				reticulations[retic] += 1
			except:
				# if reticulation not yet in dictionary, add it
				# with count 0 because we want indegree-1
				reticulations[retic] = 0
		
		retnum = 0
		for retic in reticulations:
			retnum += reticulations[retic]

		retnumlist.append(retnum)
		
	infile.close()
	return(retnumlist)

def list_files(directory, substring):
	cmd1 = f"cd {directory}"
	cmd2 = f"find *{substring}*"
	
	p1 = subprocess.run(f"{cmd1} && {cmd2}", shell=True, capture_output=True)
	filelist = p1.stdout.decode("utf8").strip().split("\n")
	
	return(filelist)

def get_network_specs(algorithm):
	outfile = open(f"{resultspath}{algorithm}_network_specs.csv", "w")
	outfile.write("Filename,Reticulation number\n")
	
	# get list of all network file names
	filelist = list_files(resultspath + algorithm, ".net")

	for filename in filelist:
		outfile.write(f"{filename},")
		retnumlist = count_reticulation_number(resultspath + algorithm + "/" + filename)
		for retnum in retnumlist:
			outfile.write(f"{retnum} ")
		outfile.write("\n")
		
	outfile.close()
	
def get_temp_specs():
	outfile = open(f"{resultspath}temporal_specs.csv", "w")
	outfile.write("Filename,Temporal distance,Reticulation number,Maximum temporal distance\n")
	
	filelist = list_files(resultspath + "temporal", ".log")
	
	for filename in filelist:
		outfile.write(f"{filename},")
		infile = open(f"{resultspath}temporal/{filename}", "r")
		
		sol = False
		pmax = ""
		for line in infile:
			if "Solution found for" in line:
				sol = True
				k = line.split()[3][2] # hybridization number of solution
				p = line.split()[5][2] # temporal distance of solution
				
			elif "maxTemporalDistance" in line:
				pmax = line.split("=")[1].strip() # temporal distance searching
				
		if sol:
			outfile.write(f"{p},{k},\n")
		else:
			outfile.write(f",,{pmax}\n")
		
		infile.close()
		
	outfile.close()

def get_nonbinary_to_binary_info(infilename):
	# output log from data transform
	# from tree-child network
	infile = open(resultspath + "nonbinary-to-binary/" + infilename, "r")
	
	for line in infile:
		if "resolved multifurcations" in line:
			line = line.split(":")
			nres = line[1].strip()
		if "output trees" in line:
			line = line.split(":")
			ntrees = line[1].strip()
	
	infile.close()
	
	return(nres, ntrees)
	
def count_multifurcations(infilename):
	# file containing non-binary trees in newick format
	# one per line
	infile = open(datapath + infilename, "r")
	nmult = 0
	nleafs = 12
	
	for line in infile:
		nleafs = len(line.split(","))
		nbrackets = line.count("(")
		nmult += (nleafs - 1 - nbrackets)
		
		#line = re.split("\(|\)", line) # split on "(" and ")"
		#for word in line:
		#	c = word.count(",") # count commas
		#	if c > 1:
		#		nmult += (c-1)
		
	infile.close()
	
	return(nmult)
	
def get_cass_specs(infilename):
	infile = open(f"{resultspath}cass/{infilename}", "r")
	
	maxlevel = ""
	error = ""
	sol = False
	for line in infile:
		if "Found" in line:
			if "trees" in line:
				ntrees = line.split()[1]
				ntax = line.split()[4]
			elif "clusters" in line:
				nclust = line.split()[1]
			elif "nontrivial" in line:
				ncomp = line.split()[1]
		elif "Searching level" in line:
			maxlevel = line.split()[2]
		elif "Finished all components!" in line:
			sol = True
		elif "Error" in line:
			error = line.strip()
	infile.close()
	
	return(ntrees, ntax, nclust, ncomp, maxlevel, sol, error)

def get_tree_specs():
	outfile = open(f"{resultspath}tree_specs.csv", "w")
	outfile.write("Version,Trees,Taxa,Clusters,Nontrivial components,ICass level,ICass maximum level searched,ICass error,Contracted edges,Resolved multifurcations,Unique trees\n")
	
	# get list of all tree names
	filelist = list_files(resultspath + "cass/", ".log")

	for filename in filelist:
		version = filename.replace("_noBS.log", "")
		outfile.write(f"{version},")
		ntrees, ntax, nclust, ncomp, maxlevel, sol, error = get_cass_specs(filename)
		if sol:
			outfile.write(f"{ntrees},{ntax},{nclust},{ncomp},{maxlevel},,{error},")
		else:
			outfile.write(f"{ntrees},{ntax},{nclust},{ncomp},,{maxlevel},{error},")
		nmult = count_multifurcations(version + ".tree")
		try:
			nres, uniquetrees = get_nonbinary_to_binary_info(version + ".log")
			outfile.write(f"{nmult},{nres},{uniquetrees}\n")
		except:
			outfile.write(f"{nmult},,\n")
		finally:
			pass
		
	outfile.close()
	
def save_png(algorithm, substring):
	filelist = list_files(resultspath + algorithm, substring)
	
	for filename in filelist:
		with open(dendropath+"temp_cmdfile.txt", "w") as cmdfile:
			cmdfile.write(f"open file='/home/rosanne/Documents/minor_research_project/{resultspath.replace('../','')}{algorithm}/{filename}';\n")
			cmdfile.write(f"exportImage file='/home/rosanne/Documents/minor_research_project/{resultspath.replace('../','')}{algorithm}/{filename.replace('{substring}','.png')}' format=png visibleOnly=false title=none replace=true;\nquit")
			
		cmd1 = "cd " + dendropath
		cmd2 = "./Dendroscope -g -c ./temp_cmdfile.txt"
		
		p1 = subprocess.call(f"{cmd1} && {cmd2}", shell=True)

#get_network_specs("treechild")
#get_network_specs("phylonet")
#get_temp_specs()
#gettree_specs()

#save_png("treechild", "seqsel10*.net")
#save_png("trilonet", ".net")
#save_png("tril2net", ".eNewick")

