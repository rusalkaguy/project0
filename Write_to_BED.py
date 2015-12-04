# Defining Accession Number, File_name, and Database

# Import adds modules to your script from the python argv
# argv = argument variable that holds the arguments you pass to your script when you run it
from sys import argv 

# Line 14 upacks argv so that it gets assigned to 1 variable you can work with
# "unpack" = take whater is in argv, unpack it, and assign it to all these variables on the left in order.
# argv = accession_number 
accession_number = argv[1]
file_name = accession_number + '.gbk'
db = 'nucleotide' 

#Write Genbank File to BED file format: https://gist.github.com/brantfaircloth/893580
from Bio import SeqIO

import pdb

def main():
    outf = open( accession_number +'.bed', 'w')
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

# Run with $ python Write_to_BED.py NC_006273