# Defining Accession Number and Database

# Import adds modules to your script from the python argv
# argv = argument variable that holds the arguments you pass to your script when you run it
from sys import argv 

# Line 14 upacks argv so that it gets assigned to 1 variable you can work with
# "unpack" = take whater is in argv, unpack it, and assign it to all these variables on the left in order.
# argv = accession_number 
accession_number = argv[1]
file_name = accession_number + '.gbk'
db = 'nucleotide'

# Grabbing genomes from Genbank
from Bio import Entrez

# Entrez sends an email reqesting the data specified below
Entrez.email = 'bheater@uab.edu'
handle=Entrez.efetch(db=db,id=accession_number,rettype='gb') 
# Accession id works, returns genbank format, looks in the 'nucleotide' database

# Store locally
local_file=open(file_name,'w') # opens and create file (W)
local_file.write(handle.read()) # write takes data and writes to file
handle.close()
local_file.close()

# run with $ python Genomes_from_Genbank.py NC_006273

'''
Create a dictionary of gene names as keys
'''
bed_gene_name_col = 3
# Open (fo = "file open") the bed file previously created.
fo = open ("NC_006273.gbk", "r+")

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
