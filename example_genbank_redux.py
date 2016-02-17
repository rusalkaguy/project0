''' 
Variables for future use writing in bed format
chrom = NC_006273 # The name of the chromosome
chromStart =  # The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chomEnd =  # The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. 
# For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
name = # Defines the name of the BED line. 
score = # A score between 0 and 1000.
strand = # Defines the strand - either '+' or '-'complement.
thickStart = # The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). 
# When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
thickEnd = # The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
itemRgb = # An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. 
blockCount = # The number of blocks (exons) in the BED line.
blockSizes = # A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
blockStart = # A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
'''

feature_type_col = 0
gene_name_col = 1
location_col = 2


# Open the file of interest
file_open = open ("example_genbank_redux.txt", "r+")
file_read = file_open.readlines()

gene_dict = {} # Creates a new dictionary
for line in file_read:
    line = line.split('#')[0]
    tab = line.split()  # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    if len(tab) == 0:
        continue
    print tab
    gene_name = tab[gene_name_col]
    feature_type = tab[feature_type_col]
    location = tab[location_col]
    if gene_name in gene_dict
'''
    # look up gene to see if it exists
    try:
		feature_list = gene_dict[gene_name] # get defition for this gene
        # If gene exists in dictionary, append new feature
        feature_list.append(feature_dict)
    except KeyError:
        # If new gene, create an array with one exon in it
        feature_list = [feature_dict]
        gene_dict[gene_name] = feature_list
'''
for key in gene_dict:
    print gene_dict[key]
