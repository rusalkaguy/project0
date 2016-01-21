'''
Create a dictionary of gene names as keys
'''
pos_col = 0
key_name_col = 1
gene_name_col = 2

# Open the file previously created with array.
file_open = open ("example_for_dictionary.py", "r+")
file_read = file_open.readlines()
# Create a new dictionary

new_dict = {} # Creates a new dictionary

for line in file_read:
    tab = line.split() # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    # When already in the dictionary
    if tab[key_name_col] in new_dict:  # Checks if the value is in the dictionary
        new_dict[tab[key_name_col]].append(str(tab[pos_col]) + ':' + str(tab[gene_name_col]))
    else:
        new_dict[tab[1]] = [str(tab[pos_col]) + ':' + str(tab[gene_name_col])]

# Print the results
for key in new_dict:
    print(str(key) + '\t' + ',' .join(new_dict[key])) 


# Close opened file
file_open.close
for item in new_dict:
    # new_dict as a dictionary has random order
    'output = str(new_dict[item]) +'    '+ str(pos_name1)'

'''
    for item in newer_dict:
        y = str(x +":" + newer_dict [item])

        print y
        # Prints all keys in random order from new_dict
'''


