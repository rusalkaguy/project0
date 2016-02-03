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
for line in file_read:
    if '#' in line:     # Change to if '#' exists in line, 
        continue        # only consider code before '#' 
    tab = line.split() # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    index = tab[index_col]
    gene_name = tab[gene_name_col]

    exon_dict = {}
    exon_dict["name"] = tab[exon_name_col]
    exon_dict["start"] = tab[exon_start_col]
    exon_dict["stop"] = tab[exon_stop_col]
    
    # look up gene to see if it exists
    exon_list = gene_dict[gene_name] # Sets an empty dictionary
    if gene_name not in exon_list: # KeyError: 'A'
        # gene exists, just append new exon
        exon_list.append(exon_dict) 
    else:
        # new gene - create an array with one exon it it
        gene_dict[gene_name] = exon_dict
    
# Shift print formatting from prior file

print "#--------- gene_dict ------------"
for key in gene_dict:
    # create string for exon list
    gene_ex_list = gene_dict[key]
    ex_output_array = []
    for ex in gene_ex_list:
        # format exon (ex) and add to output array
        formatted_exon_string = ex[1] + ":" + ex[2] + "-" + ex[3]
        ex_output_array.append(formatted_exon_string)
        #print "%", key, ":" , formatted_exon_string
        #if ex not in ex_output_array:
        #    ex_output_array.append(ex[1])
        #    ex_output_array.append(ex[2:4])
        #    continue
        #if ex in ex_output_array:
        #    ex_output_array.append(ex[2:4])
        #    continue
    #print ex_output_array
    print key+'\t'+','.join(ex_output_array)

# Print the results
for key, value in gene_dict.iteritems():
    print(key, value) 

print gene_dict
# Understand flow and formatting: draw a picture.

# Close opened file
file_open.close