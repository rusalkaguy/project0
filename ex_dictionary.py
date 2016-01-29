'''
Create a dictionary of gene names as keys
'''
index_col = 0
gene_name_col = 1
exon_name_col = 2
exon_start_col = 3
exon_stop_col = 4


# Open the file previously created with array.
file_open = open ("example_for_dictionary.txt", "r+")
file_read = file_open.readlines()
# Create a new dictionary

gene_dict = {} # Creates a new dictionary
exon_dict = {}
for line in file_read:
    if '#' in line:
        continue
    tab = line.split() # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    index = tab[index_col]
    gene_name = tab[gene_name_col]
    exon_name = tab[exon_name_col]
    exon_start = tab[exon_start_col]
    exon_stop = tab[exon_stop_col]
    exon_list = []
    # When the gene is not present in the dictionary
    if gene_name not in gene_dict:  # Checks if the gene is not in the dictionary
        gene_dict[gene_name] = exon_list.append(exon_name)
    else:
        gene_dict[gene_name] = exon_list.append(exon_name)
    for exon_name in exon_list:
        exon_dict [exon_name] = [index, gene_name, exon_name,exon_start,exon_stop]


# Print the results
for key, value in gene_dict.iteritems():
    print(key, value) 
for key, value in exon_dict.iteritems():
    print(key, value)
# Close opened file
file_open.close