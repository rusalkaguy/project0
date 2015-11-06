#Write Genbank File to GFF file format: http://biopython.org/wiki/GFF_Parsing
from BCBio import GFF
from Bio import SeqIO
 
in_file = "your_file.gb"
out_file = "your_file.gff"
in_handle = open(in_file)
out_handle = open(out_file, "w")
 
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
 
in_handle.close()
out_handle.close()
