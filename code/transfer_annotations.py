path = "../data/"
mesquitefile = "MSA_annotation_deleted_regions.txt"
raxmlpath = "../ubuntu_software/standard-RAxML-master/"
originalmsafile = "Corona_original_big.mafft.nex"

def make_delregion_list(infilename):
	# from mesquite file
	with open(infilename, "r") as infile:
		charpart = infile.read().split('\n')[4]
	charpart = charpart.split("=")[1].split(",")
	delregion_list = []

	for delregion in charpart:
		delregion = delregion.split("-")
		start = int(delregion[0].strip())
		end = int(delregion[1].strip()) # note that end is t/m so we do +1
		length = end + 1 - start
		delregion_list.append([start,end+1,length])

	delregion_list.sort()

	return(delregion_list)
	
def get_gap_locations(originalmsafilename):
	# get SARS-CoV-2 gap locations from the original MSA file
	infile = open(originalmsafilename)
	gaplocations = []
	for line in infile:
		if line.strip().startswith("'NC_045512_Wuhan_Hu_1'"):
			seq = line.split()[1].strip()
			for c in range(0, len(seq)):
				if seq[c] == "-":
					if seq[c-1] != "-":
						startgap = c+1 # index starts at 1 so +1
					if seq[c+1] != "-":
						endgap = c+2 # we want 'tot' not 't/m' so +2
						gaplocations.append([startgap, endgap])
			break
					
	infile.close()
	
	return(gaplocations)

def update_start_end_gaps(gap, start, end):
	# transfer SARS-CoV-2 gene annotation to mesquite alignment
	if gap[0] <= start:
		start += (gap[1]-gap[0])
		end += (gap[1]-gap[0])
	elif gap[0] < end:
		end += (gap[1]-gap[0])
			
	return(start, end)

def move_start_end_deletions(delregion, start, end):
	# transfer mesquite alignment to adjusted MSA
	mvstart = 0
	mvend = 0
	
	if delregion[0] < start:
		if delregion[1] < start:
			# entire delregion before the gene		
			mvstart = (delregion[1]-delregion[0])
		elif delregion[1] < end:
			# delregion starts before gene and ends within gene
			print(f"delregion starts before and ends within gene")
			mvstart = (start-delregion[0])
		else:
			# delregion spans entire gene
			print(f"Entire gene deleted")
		mvend = (delregion[1]-delregion[0])
			
	elif delregion[0] < end:
		if delregion[1] < end:
			# entire delregion within gene
			mvend = (delregion[1]-delregion[0])
		else:
			# delregion starts within gene and ends after gene
			print(f"delregion starts within and ends after gene")
			mvend = (end-delregion[0])
	
	return(mvstart, mvend)

def transfer_annotations(genefilename, mesquitefilename, originalmsafilename, outfilename):
	with open(genefilename, "r") as infile:
		gene_list = infile.read().strip().split("\n")
		
	delregion_list = make_delregion_list(path+mesquitefilename)
		
	gaplocations = get_gap_locations(path+originalmsafilename)

	outfile = open(path+outfilename, "w")
	# SARS-CoV-2 has a 38 base 'gap' at the beginning of the 
	# mesquite and MSA files so all genes have to be shifted 
	# +38 bases to align with the delregions
	shift = 38

	for gene in gene_list:
		gene = gene.split()
		name = gene[0]
		
		start = int(gene[1]) + shift
		end = int(gene[2]) + 1 + shift # +1 because we want 'tot' not 't/m'
		
		# transfer SARS-CoV-2 gene annotation to mesquite alignment
		for gap in gaplocations:
			start, end = update_start_end_gaps(gap, start, end)
			
		# transfer mesquite alignment annotation to adjusted MSA
		newstart = start
		newend = end
		for delregion in delregion_list:
			mvstart, mvend = move_start_end_deletions(delregion, start, end)
			newstart -= mvstart
			newend -= mvend

		outfile.write(f"{name} {newstart} {newend}\n")
	
	outfile.close()
	
def check_genes(checkfilename, finalmsafile, version):
	with open(checkfilename, "r") as checkfile:
		gene_list = checkfile.read().strip().split("\n")
		
	for gene in gene_list:
		name = gene.split()[0]
		infile = open(f"{path}{finalmsafile.replace('.ph','')}_{version}_dir/{finalmsafile}.{version}.{name}.txt", "r")
		for line in infile:
			if line.strip().startswith("MN908947_Wuhan_Hu_1"):
				seq = line.strip().split()[1]
				if seq.startswith("ATG"):
					print(f"{name} starts with start codon ATG")
				else:
					print(f"Check this gene: {name} starts with {seq[0:40]}")
		infile.close()
		
def info_on_gene_annotations(annotationfilename, outfilename):
	annotationfile = open(path+annotationfilename, "r")
	outfile = open(path+outfilename, "r")
	
	for line in annotationfile:
		line = line.strip().split()
		print(f"Gene {line[0]} has length {int(line[2])-int(line[1])} originally in SARS-CoV2")
		line = outfile.readline().strip().split()
		print(f"Gene {line[0]} has length {int(line[2])-int(line[1])-1} in the MSA file") # minus one because 't/m' changed to 'tot'

annotationfile = "gene_annotations_all_sarscov2.txt"
outfile = "gene_annotations_all_MSA.txt"

#transfer_annotations(path+annotationfile, path+mesquitefile, path+originalmsafile, path+outfile)

#info_on_gene_annotations(annotationfile, outfile)

outfile = "gene_annotations_selection_MSA.txt"
finalmsafile = "Corona_edited_seqsel7.ph"
version = "genes_seqsel7"
check_genes(path+outfile, finalmsafile, version)
