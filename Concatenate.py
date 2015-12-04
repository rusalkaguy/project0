# This program cooncatenates Genomes_from_Genbank.py with Write_to_BED.py and Write_to_GFF.py into one file.

# Defining Accession Number and Database

# Import adds modules to your script from the python argv
# argv = argument variable that holds the arguments you pass to your script when you run it
from sys import argv 

# Line 14 upacks argv so that it gets assigned to 1 variable you can work with
# "unpack" = take whatever is in argv, and assign it to all these variables on the left in order.
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

#Write Genbank File to BED file format: https://gist.github.com/brantfaircloth/893580
from Bio import SeqIO

import pdb

def main():
    outf = open( accession_number +'.bed', 'w')
    header = accession_number
    outf.write(header)
    for record in SeqIO.parse(open(file_name, "rU"), "genbank") :
        for feature in record.features:
            if feature.type == 'gene':
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = feature.qualifiers['gene'][0]
                except:
                    # some features only have a locus tag
                    name = feature.qualifiers['locus_tag'][0]
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = accession_number + "\t{0}\t{1}\t{2}\t1000\t{3}\t{0}\t{1}\t65,105,225\n".format(start, stop, name, strand)
                outf.write(bed_line)
    outf.close()


if __name__ == '__main__':
    main()

# Write Genbank File to GFF file format: http://biopython.org/wiki/GFF_Parsing
from BCBio import GFF
from Bio import SeqIO
 
in_file = file_name
out_file = accession_number + '.gff'
in_handle = open(in_file)
out_handle = open(out_file, "w")
 
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
 
in_handle.close()
out_handle.close()

# Concatenate files written in BED and GFF into one file.
filenames = ['accession_number.bed', 'accession_number.gff', ...]
with open('path/to/output/file', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)