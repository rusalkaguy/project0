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