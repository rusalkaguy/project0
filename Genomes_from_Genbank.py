#Defining Accession Number and Database
accession_number= 'U00096'
file_name = accession_number + '.download.gbk'
db='nucleotide' #To determine database from assession prefixes:http://www.ncbi.nlm.nih.gov/Sequin/acc.html

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
#clean up resources