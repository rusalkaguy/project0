'''
Create a dictionary of gene names as keys
'''
bed_gene_name_col = 3
# Open (fo = "file open") the bed file previously created.
fo = open ("gbout.bed", "r+")

# Create a new dictionary og gene names
gene_dict ={} # Key is the gene name 

for line in fo:
    tab = line.split() # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    gene_name = tab[bed_gene_name_col]
    gene_dict[gene_name] = tab
    # Brackets denote the [key] is the gene name in collumn 3
    # The '' assigns a null value to the keys in the dictionary.
    # newdict[tab[3]] could =  tab or line if all data is valuable.
    # Tab is more helpful becuase tab is already parsed.

# Close opened file
fo.close

for item in gene_dict:
    # Gene_dict as a dictionary has random order
    print item, gene_dict[item]
    # Prints all keys (gene names) in random order from gene_dict
