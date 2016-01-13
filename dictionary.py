'''
Create a dictionary of gene names as keys
'''
pos_col = 0
key_name_col = 1
pos_name_col = 2

# Open (fo = "file open") the file previously created with array.
fo = open ("example_for_dictionary.py", "r+")
fr = fo.readlines()
# Create a new dictionary
new_dict = {} # Key is the letter in the first collumn. 

for line in fr:
    tab = line.split() # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    key_name = tab[key_name_col]
    if key_name in new_dict:
        
    else new_dict[key_name] = tab [pos_col]
    # Brackets denote the [key] is the key_name
    # new_dict[key_name] = tab [pos_col] sets the key name to be the position
    # Tab is more helpful becuase tab is already parsed.

    # Create a nested dictionary
    newer_dict = {}

    for line in fr:
        tab = line.split()
        pos_name = tab[pos_col]
        newer_dict[pos_name] = tab[pos_name_col]

# Close opened file
fo.close

for item in new_dict:
    # new_dict as a dictionary has random order
    x = str(item +"    "+ new_dict[item])


    for item in newer_dict:
        y = str(x +":" + newer_dict [item])

        print y
        # Prints all keys in random order from new_dict

