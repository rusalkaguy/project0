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

# Creates an array for consisting of each line of the file
file_read = file_open.readlines()
# file_read = ['line1','line2',...]

gene_exon_array=[]

for line in file_read:
    # Ignore first line with comment
    if '#' in line:
        continue
    # Split the line in the file by tab

    # Create an array of items within each line
    row = line.split()
    gene_exon_array.append(row)
    #                       line1          line 2
    # gene_exon_array = [['1','A','qwe','22','3'], ['2','B','wse','1','4'],...]

# Create the array for gene names
gene_name_array = []

for line in gene_exon_array:
    if line [1] not in gene_name_array:
        gene_name_array.append(line[1])
        # gene_name_array = ['A','B','C','D','Z']

# Create the dictionary for exons
exon_dict = {}
# Create the dictionary for gene names
gene_dict = {}

for gene in gene_name_array:
    exon = []
    for line in gene_exon_array:
        if line [1] == gene:
            exon.append(line[0:1]+line[2:6])
    # Set dictionary key to be gene correlated to exon information
    gene_dict[gene]=exon


    for item in exon:
        if item[0] not in exon_dict:
            exon_dict[item[0]]=item[1]


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

print "#--------- gene_dict ------------"
pp.pprint(gene_dict)

# Close opened file
file_open.close

'''
exon_list = [] # Initializes an empty array (not a list)
# Create a new dictionary
gene_dict = {} # Creates a new dictionary
for line in file_read:
    if '#' in line:
        continue
    index = tab[index_col]
    gene_name = tab[gene_name_col]
    exon_name = tab[exon_name_col]
    exon_start = tab[exon_start_col]
    exon_stop = tab[exon_stop_col]
    exon_dict = {}
    ### put name/value pairs into exon_dict

    exon_dict[exon_name] = [index, gene_name, exon_name,exon_start,exon_stop]
    # When the gene is not present in the dictionary
    if gene_name not in gene_dict:  # Checks if the gene is not in the dictionary
        
        exon_list.append(exon_dict)
        gene_dict[gene_name] = exon_list
    else:
        exon_list=gene_dict[gene_name]
        exon_list.append(exon_dict)
    # Unneccessary:
    # for exon_name in exon_list:
    #    exon_dict [exon_name] = [index, gene_name, exon_name,exon_start,exon_stop]


# Print the results
for key, value in gene_dict.iteritems():
    print(key, value) 
    exon_list = value;
    for item in exon_list:
        exon_dict[exon_name]=exon_list
        for key, value in exon_dict.iteritems():
            print(key, value)
for i in gene_dict:
    print i 
    raw_input()

'''

