#!/usr/local/bin/python
# coding: latin-1
import os, sys
'''This script should 
	1. read in a bed file, 
	2. break up the 4th column for gene name into gene_name and genome_accession, 
	3. reformat the bed file with genome_accession as an additional, searchable identifier, and 
	4. create file and directories for the track hub'''
# Use argv to obtain input of accession # from command line
from sys import argv 
from os import path
import errno
import os
import shutil

bed_file_name = argv[1] # Upacks argv-> assigned to 1 variable you can work with
abrev = argv[2] # abrev is the the abreviated name of the chrom.sizes file and subdirectory, which should be SR10-01
genome = abrev
track_hub_directory = 'hcmv_pub'

# Read in bed file
bed_file=open(bed_file_name,'r')
line_list =bed_file.readlines()
col_list = line_list[1].split('\t')
# accession_number = col_list[0]
# genome = accession_number
#print line_list
num_lines = len(line_list)
line = 1
processed_bed_file_name= 'processed_'+bed_file_name
processed_bed_file = open(processed_bed_file_name,'w')
output = []
while line<num_lines:
	col_list = line_list[line].split('\t')
	#print col_list
	#print 'length of col_list:', len(col_list)
	#print line
	chrom = col_list[0]
	chromStart = col_list[1]
	chromEnd = col_list[2]
	gene_name = col_list[3]
	# gene_name_and_genome = col_list[3].split('.')
	# gene_name = gene_name_and_genome[0]
	# genome_name = gene_name_and_genome[1]
	#print gene_name +','+ genome_accession
	score = col_list[4]
	strand = col_list[5]
	thickStart = col_list[6]
	thickEnd = col_list[7]
	itemRgb = col_list[8]
	blockCount = col_list[9]
	blockSizes = col_list[10]
	blockStarts = col_list[11]
	new_line = '\t'.join([chrom,chromStart,chromEnd,gene_name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts])
	output.append(new_line)
	#print output
	#raw_input()
	line = line+1
#print output
bed_file.close
for line in output:
	processed_bed_file.write(line)
processed_bed_file.close



def mkdir_p(abrev):
	try:
		os.makedirs(abrev)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(abrev):
			pass
		else:
			raise
	print "New directory, "+abrev+", created in the current directory."

def sort_bed_file(bed_file_name, genome):
	old_filename = bed_file_name
	bed_file = open(bed_file_name)
	new_filename = genome+'sorted.bed'
	from subprocess import call
	cmd = ["sort","-k1,1","-k2,2n", old_filename,'-o', new_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)

def copy_file(filename,track_hub_directory,abrev):
	subdir2 = abrev
	subdir1 = track_hub_directory
	dest_filepath = os.path.join(subdir1, subdir2, filename)
	try:
		shutil.copy(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir1+'\\'+subdir2
	except IOError:
		print "Wrong path provided because the bed file does not exist."

def mk_chrom_sizes_file(genome,track_hub_directory,abrev):
	# 1. If known database, $ fetchChromSizes <db> > <db>.chrom.sizes
	# 2. If new database, determine the number of nucleotides in the chromosomes
	subdir = track_hub_directory
	filename = abrev+".chrom.sizes"
	dest_filepath = os.path.join(subdir, filename)
	file_contents = genome+"\t235646"
	chrom_sizes_file = open(dest_filepath,'w')
	chrom_sizes_file.write(file_contents)
	print 'wrote '+filename+' to '+dest_filepath
	chrom_sizes_file.close
	chrom_sizes_file = open(filename,'w')
	chrom_sizes_file.write(file_contents)
	chrom_sizes_file.close
	print 'wrote '+filename+' to current directory'

def bedToBigBed(genome,abrev):
	from subprocess import call
	sorted_bed_file = genome+'sorted.bed'
	output_filename = genome+'.bb'
	chrom_sizes_file = abrev+'.chrom.sizes'
	# bedToBigBed -extraIndex=name -type=bed12 fileinsorted.bed chrom.sizes outputfile.bb
	index = '-extraIndex=name'
	bed_type = '-type=bed12'
	cmd = ["bedToBigBed",index,bed_type, sorted_bed_file, chrom_sizes_file, output_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)

def mk_genomes_file(genome,track_hub_directory):
	filename = 'genomes.txt'
	# genomes_str = 'genome '+abrev+'\ntrackDb '+abrev+'/trackDb.txt\ntwoBitPath '+abrev+'/NC_006273v2.2bit\ngroups '+abrev+'/groups.txt\ndescription '+genome+'\norganism HH5 strain Merlin\ndefaultPos NC_006273v2:1-235646\norderKey 100\nscientificName Human herpesvirus 5\nhtmlPath description.html'
	# genomes_file = open(filename,'w')
	# genomes_file.write(genomes_str)
	# genomes_file.close
	line1 = 'genome '+abrev
	line2 = 'trackDb '+abrev+'/trackDb.txt'
	line3 = 'twoBitPath '+abrev+'/'+genome+'.2bit'
	line4 = 'groups '+abrev+'/groups.txt'
	line5 = 'description '+genome
	line6 = 'organism HH5 strain Merlin'
	line7 = 'defaultPos '+genome+':1-235646'
	line8 = 'orderKey 100'
	line9 = 'scientificName Human herpesvirus 5'
	line10 = 'htmlPath description.html'
	genomes_list_of_lines = [line1,line2,line3,line4,line5,line6,line7,line8,line9,line10]
	genomes_file_open = open(filename,'w')
	for line in genomes_list_of_lines:
		genomes_file_open.writelines(line+'\n')
	genomes_file_open.close()

	# Save file to subdir
	subdir = track_hub_directory
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

def change_dir(track_hub_directory):
	try:
		os.chdir(track_hub_directory)
	except WindowsError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(track_hub_directory):
			pass
		else:
			raise

def mk_groups_file(genome,abrev):
	filename = 'groups.txt'
	groups_str = '''name gene\nlabel Genes and Gene Predictions\npriority 10\ndefaultIsClosed 1'''
	groups_file = open(filename,'w')
	groups_file.write(groups_str)
	groups_file.close

def mktrackDb_file(genome,abrev):
	filename = "trackDb.txt"
	trackDb_str = "track "+abrev+"_refseq_mrna\nbigDataUrl "+genome+".bb\nshortLabel RefSeq Transcripts\nlongLabel RefSeq transcripts(mRNA)\ncolorByStrand 150,100,30 230,170,40\ncolor 150,100,30\naltColor 230,170,40\ntype bigBed 12\ngroup genes\nsearchIndex name"
	trackDb_file = open(filename,'w')
	trackDb_file.write(trackDb_str)
	trackDb_file.close
'''
def fasta_file(accession_number, abrev):
	# if accession number has a period, indicating version number, replace with v.
	if '.' in accession_number:
		genome='v'.join(accession_number.split('.'))
	else:
		genome = str(accession_number)
	file_name = genome + '.fna'
	fasta_file = open(file_name, 'w')
	from Bio import SeqIO
	fasta_contents = SeqIO.convert(accession_number+".gbk", "genbank", genome+".fna", "fasta")
	# add code here to edit first line of fasta file_content by changing accession_number to path_file
	# parsing code
	fasta_file.close
	fasta_file=open(file_name,'r')
	fasta_list =fasta_file.readlines()
	#print fasta_list[0]
	words = fasta_list[0].split(' ')
	if '.' in words[0]:
		words[0]='>'+genome
		fasta_list[0]=' '.join(words)
		fasta_text= ''.join(fasta_list)
	else:
		fasta_text=fasta_list
	fasta_file.close
	fasta_file = open(file_name, 'w')
	fasta_file.write(fasta_text)
	fasta_file.close
'''
def fasta_to_2bit(fasta_file):
	from subprocess import call
	#fasta_file = genome+'.fna'
	output_2bit_filename = genome+'.2bit'
	cmd = ["faToTwoBit",fasta_file, output_2bit_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)

if __name__ == '__main__':
	mkdir_p(track_hub_directory)
	sort_bed_file(processed_bed_file_name,genome)
	copy_file(genome+"sorted.bed",track_hub_directory,abrev)
	mk_chrom_sizes_file(genome,track_hub_directory,abrev)
	bedToBigBed(genome,abrev)
	copy_file(genome+'.bb',track_hub_directory,abrev)

	mk_hub_txt_file(genome,track_hub_directory)
	mk_descriptionUrl_file(genome,track_hub_directory)
	mk_genomes_file(genome,track_hub_directory)
	#copy_file(accession_number+'.gbk',track_hub_directory,abrev)
	# dont have or need a genbank file because its new data
	fasta_to_2bit(genome,abrev)
	# copy 2 bit file to subdir
	copy_file(genome+'.2bit',track_hub_directory,abrev)

	change_dir(track_hub_directory)
	mkdir_p(abrev)
	change_dir(abrev)
	# check locations of file created in functions below
	mk_groups_file(genome,abrev)
	mktrackDb_file(genome,abrev)
	#fasta_file(accession_number,abrev)

	

# to run
# python bed_to_trackhub.py all_geno.SR10-01.vicuna.contig.score100.best1.bed SR10-01