'''
Create a dictionary of gene names as keys
'''
index_col = 0
gene_name_col = 1
exon_name_col = 2
exon_start_col = 3
exon_stop_col = 4

# For debugging output
import pprint
pp = pprint.PrettyPrinter(indent=4)

# Open the file.
file_open = open ("example_for_dictionary.txt", "r+")
file_read = file_open.readlines()
# Create a new dictionary

gene_dict = {} # Creates a new dictionary

line = file_read.splitlines()
gene_map = map(lambda line: line.split()[exon_name_col], line)
'''
for line in file_read:
    if '#' in line:     # Change to if '#' exists in line, 
        continue        # only consider code before '#' 
    tab = line.split()  # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    index = tab[index_col] # Not neccessary anymore, but may be helpful later with g
    gene_name = list(tab[gene_name_col]) # Set gene_name to be the item in the gene_name column.
    exon_name = tab[exon_name_col]
    exon_start = tab[exon_start_col]
    exon_stop = tab[exon_stop_col]
    gene_map = map(None,gene_name)
    exon_map = map(None,exon_name)
    exon_start_map = map(None,exon_start)
    exon_stop_map = map(None,exon_stop)
    gene_exon_map = "\t".join(gene_map)
    print gene_exon_map


    exon_dict = {}
    exon_dict["name"] = tab[exon_name_col]
    exon_dict["start"] = tab[exon_start_col]
    exon_dict["stop"] = tab[exon_stop_col]
    
    # look up gene to see if it exists
    try: 
        exon_list = gene_dict[gene_name] # get defition for this gene
        # gene exists, just append new exon
        exon_list.append(exon_dict) 
    except KeyError:
        # if a new gene, create an array with one exon in it
        exon_list = [exon_dict]
        # set gene_dict with key of gene_name 
        # to be exon_list, the array with exon_dict,
        gene_dict[gene_name] = exon_list


print "#--------- gene_dict ------------"
pp.pprint(gene_dict)
print "#--------- formatted gene_dict ------------"
for key in gene_dict:
    # create string for exon list
    gene_ex_list = gene_dict[key]
    ex_output_array = []
    for exon_dict in gene_ex_list:
        # format exon (exon_dict) and add to output array
        formatted_exon_string = exon_dict['name'] + ":" + exon_dict['start'] + "-" + exon_dict['stop']
        ex_output_array.append(formatted_exon_string)
    print key+'\t'+','.join(ex_output_array)
'''

# Close opened file
file_open.close