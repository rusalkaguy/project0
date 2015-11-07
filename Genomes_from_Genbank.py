
#Defining Accession Number and Database
accession_number= 'U00096'
file_name = accession_number + '.download.gbk'
db='nucleotide' #To determine database from assession prefixes:http://www.ncbi.nlm.nih.gov/Sequin/acc.html

from sys import argv # Import adds modules to your script from the python argv; argv = argument variable that holds the arguments you pass to your script when you run it
argv= accession_number  #Line 8 upacks argv so that it gets assigned to 1 variable you can work with; "unpack" = take whater is in argv, unpack it, and assign it to all these variables on the left in order.

#Grabbing genomes from Genbank
from Bio import Entrez
#Replace with your real email 
Entrez.email = 'bheater@uab.edu'
handle=Entrez.efetch(db=db,id=accession_number,rettype='gb') # Accession id works, returns genbank format, looks in the 'nucleotide' database
#store locally
local_file=open(file_name,'w') #opens and create file (W)
local_file.write(handle.read())#write takes data and writes to file
handle.close()
local_file.close()

