'''
Create a dictionary of gene names as keys
'''
index_col = 0
gene_name_col = 1
exon_name_col = 2
exon_start_col = 3
exon_stop_col = 4

# for debugging output
import pprint
pp = pprint.PrettyPrinter(indent=4)

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
    gene_name = tab[gene_name_col]

    # Create exon dictionary
    exon_dict = {}
    # Set key-value pairs in exon dictionary
    exon_dict["name"] = tab[exon_name_col]
    exon_dict["start"] = tab[exon_start_col]
    exon_dict["stop"] = tab[exon_stop_col]
    
    # look up gene to see if it exists
    try: 
        exon_list = gene_dict[gene_name] # get defition for this gene
        # If gene exists in dictionary, append new exon
        exon_list.append(exon_dict) 
    except KeyError:
        # If new gene, create an array with one exon in it
        exon_list = [exon_dict]
        gene_dict[gene_name] = exon_list
    
# Shift print formatting from prior file
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
        #print "%", key, ":" , formatted_exon_string
        #if exon_dict not in ex_output_array:
        #    ex_output_array.append(exon_dict[1])
        #    ex_output_array.append(exon_dict[2:4])
        #    continue
        #if exon_dict in ex_output_array:
        #    ex_output_array.append(exon_dict[2:4])
        #    continue
    #print ex_output_array
    print key+'\t'+','.join(ex_output_array)

# Replace inner for loop above and ex_ouput_array through line 55 with map and function 
print "#--------- formatted gene_dict with map and function ------------"
def format_exon(x): # x is input argument
    return x["name"]+':'+str(x["start"])+'-'+str(x["stop"])
for key in gene_dict:
    exon_list = gene_dict[key]
    print key+'\t'+",".join(map(format_exon,exon_list))
    #for exon in exon_list:
    #    print key+'\t'+ format_exon(exon)

'''
for key in gene_dict:
    # create string for exon list
    gene_ex_list = gene_dict[key]
    def format_exon(exon_list):
        return exon_list["name"]+':'+str(exon_list["start"])+'-'+str(exon_list["stop"])
    # print gene_dict
    print format_exon(gene_ex_list)

# "name:start-stop"
print format_exon(gene_ex_list) 

# print exons lists
for xl in [gene_ex_list]:
    print xl
    print ",".join(map(format_exon,xl))
'''
# Close opened file
file_open.close