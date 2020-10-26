import subprocess

# Set path names
datapath = "../data/"
phylonetpath = "../ubuntu_software/phylonet/"
casspath = "../ubuntu_software/cass/"
dendropath = "../ubuntu_software/dendroscope/"
treechildpath = "../ubuntu_software/tree_child_code-master/"
resultspath = "/home/rosanne/Documents/minor_research_project/results/"
raxmlpath = "../ubuntu_software/standard-RAxML-master/"

def exclude_sequences(msa_fname, excl_seqs_fname):
	"""
	Works only for Phylip format.
	"""
	infile = open(datapath + msa_fname, "r")
	outfile = open(datapath + msa_fname.split(".")[0] + "_" + excl_seqs_fname.replace(".txt",".ph"), "w")
	
	with open(datapath + excl_seqs_fname, "r") as excl_seqs_file:
		sequences = excl_seqs_file.read().strip().split("\n")
	
	firstline = infile.readline().split()
	ntax = int(firstline[0].strip())
	nchar = int(firstline[1].strip())
	outfile.write(f"{ntax-len(sequences)} {nchar}\n")
	
	for line in infile:
		seq = line.split()[0].strip()
		if seq in sequences:
			continue
		else:
			outfile.write(line)
			
	infile.close()
	outfile.close()

def get_nchar_phylip(msa_fname):
	infile = open(datapath+msa_fname,"r")
	line = infile.readline().strip().split()
	nchar = int(line[1])
	infile.close()
	return(nchar)

def get_breakpoints(gene_fname):
	breakpoints = []
	infile = open(datapath + gene_fname,"r")
	for line in infile:
		line = line.strip().split()
		breakpoints.append([line[0],int(line[1]),int(line[2])])
		
	infile.close()
	return(breakpoints)

def create_excludefiles(version, msa_fname, breakpoints):
	nchar = get_nchar_phylip(msa_fname)
	# with named breakpoints
	for bp in breakpoints:
		name = bp[0]
		start = bp[1]
		end = bp[2]
		file = open(raxmlpath+version+f".{name}.txt", "w")
		if start == 1:
			file.write(f"{end}-{nchar}")
		elif end == nchar:
			file.write(f"1-{start-1}")
		else:
			file.write(f"1-{start-1} {end}-{nchar}")
		file.close()

	print("Exclude files for RAxML have been written.")
	
def run_raxml(version, msa_fname, breakpoints):
	# requires Phylip formatted MSA
	# with named breakpoints

	for bp in breakpoints:
		name = bp[0]
		cmd1 = f"cd {raxmlpath}"
		cmd2 = f"./raxmlHPC-AVX -m GTRCAT -E {version}.{name}.txt -n {version}.{name} -s ../{datapath}{msa_fname}"
		cmd3 = f"./raxmlHPC-PTHREADS-AVX -f a -m GTRCAT -s ../{datapath}{msa_fname}.{version}.{name}.txt -n {version}.{name} -p 123 -T 8 -x 123 -# 100 -w {resultspath}raxml"
		try:
			p1 = subprocess.run(cmd1+" && "+cmd2+" && "+cmd3, shell=True)
		except:
			print(f"An error has occured while trying to build a tree for {name}")
		finally:
			pass
			
	# Move the separate gene MSA's to another directory
	cmd1 = f"cd /home/rosanne/Documents/minor_research_project/data"
	cmd2 = f"mkdir {msa_fname.replace('.ph','')}_{version}_dir"
	cmd3 = f"mv ./{msa_fname}.{version}*.txt* ./{msa_fname.replace('.ph','')}_{version}_dir"
	p2 = subprocess.run(cmd1+" && "+cmd2, shell=True)
	p3 = subprocess.run(cmd1+" && "+cmd3, shell=True)
		
def cat_raxml_trees(version, name_fname):
	"""
	Concatenates RAxML trees without and with bootstrap support values
	in the order of the names listed in the namefile.
	"""
	namefile = open(datapath + name_fname, "r")
	outfile1 = open(datapath + version + ".tree", "w")
	outfile2 = open(datapath + version + "_BS.tree", "w")
	
	for line in namefile:
		name = line.strip().split()[0]
		
		try:
			with open(f"{resultspath}raxml/RAxML_bestTree.{version}.{name}","r") as infile1:
				outfile1.write(infile1.read())
			with open(f"{resultspath}raxml/RAxML_bipartitions.{version}.{name}","r") as infile2:
				outfile2.write(infile2.read())
				
		except:
			print(f"Could not find trees for {name}")
			
		finally:
			pass
	
	outfile1.close()
	outfile2.close()
	namefile.close()

"""	
def get_identical_seqs(info_fname):
	# from RAxML info file
	infofile = open(info_fname, "r")
	for line in infofile:
		if ("IMPORTANT WARNING: Sequences" in line) and ("exactly identical" in line):
			seq1 = line.split()[3]
			seq2 = line.split()[5]
			print(f"Sequences {seq1} and {seq2} are identical")
	infofile.close()
		
nchar = get_nchar_phylip(msa_fname)
breakpoints = get_breakpoints(gene_fname)
create_excludefiles(breakpoints)
run_raxml(breakpoints)
cat_raxml_trees(gene_fname)

genefile = open(datapath+gene_fname, "r")
for line in genefile:
	name = line.strip().split()[0]
	print(f"Checking info for gene {name}")
	get_identical_seqs(f"{resultspath}raxml/RAxML_info.{version}.{name}")
	
genefile.close()
"""

