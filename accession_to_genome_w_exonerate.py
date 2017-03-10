#!/usr/local/bin/python
# coding: latin-1
# github test
import os, sys
'''
This script should 
	1. Use argv to obtain input of accession # from command line, 
	2. Create directory, named by accession #, to place relevant files
	3. Run the genbank_to_bed12.py script,
		a. Pull genbank file based on accession #,
		b. parse genbank file by line, and 
		c. output .bed file using genbank_to_bed12.py; 
	4. Save output .bed file in directory and on cheaha server.
	5. Automate conversion from bed file to bigbed file on cheaha
		a. Create chrom.sizes file
			1. If known database, $ fetchChromSizes <db> > <db>.chrom.sizes
			2. If new database, determine the number of nucleotides in the chromosomes
		b. Sort the bed file: $ sort -k1,1 -k2,2n unsorted.bed > sorted.bed
		c. Convert bed to bigbed
			1. $ module load ngs-ccts/ucsc_kent/2014-03-05
			2. $ bedToBigBed in.bed hg19.chrom.sizes out.bb
		d. Save bigbed file in directory
	6. Create the hub.txt file and save in directory on server
		hub hub_name 
		shortLabel hub_short_label
		longLabel hub_long_label
		genomesFile genomes_filelist
		email email_address
		descriptionUrl descriptionUrl
	7. Create the genomes.txt file  
		a. Note: trackDb - the relative path of the trackDb file for the assembly designated by the genome tag.
			the trackDb file should be located in a subdirectory of the hub directory (preferred)
			the trackDb tag may also specify a complete URL
		b. save on directory(server)
	8. Create the trackDb.txt files
		a. Minimum Requirements: 
			track track_name
			bigDataUrl track_data_URL
			shortLabel short_label
			longLabel long_label
			type track_type
		b. Track: symbolic name of track 
			must be unique
			must be the first entry in the trackDb.txt file
		c. BigDataUrl: the file name, path, or Web location of the track's data file
			Can be URL
			path relative to the trackDb.txt file is preferable because changes to the server name 
			won’t affect relative path, but will render the URL ineffective
		d. shortLabel: short name for track appears in left column of the genome browswer
			Max length 17
		e. longLabel
		“type bigBed #” where # is the number of features or columns in the bed file
	9. Create track description files for each track and save in subdirectory
		Purpose: provide information about research to build credibility of data 
		and help people decide whether or not to use data
'''


# Create directory, named by accession #, to place relevant files

# Use argv to obtain input of accession # from command line
from sys import argv 
from os import path
import errno
import os
import shutil

accession_number = argv[1] # Upacks argv-> assigned to 1 variable you can work with
abrev = argv[2]

track_hub_directory = 'hcmv_pub'

# if accession number has a period, indicating version number, replace with v.
if '.' in accession_number:
	genome='v'.join(accession_number.split('.'))
else:
	genome = str(accession_number)
#print genome

def mkdir_p(abrev):
	try:
		os.makedirs(abrev)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(abrev):
			pass
		else:
			raise
	print "New directory, "+abrev+", created in the current directory."

def change_dir(track_hub_directory):
	try:
		os.chdir(track_hub_directory)
	except WindowsError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(track_hub_directory):
			pass
		else:
			raise

def mksubdir(abrev):
	try:
		os.makedirs(abrev)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(abrev):
			pass
		else:
			raise
	print "New directory, "+abrev+", created in the current directory."
	
def copy_file(filename,track_hub_directory,abrev):
	subdir2 = abrev
	subdir1 = track_hub_directory
	dest_filepath = os.path.join(subdir1, subdir2, filename)
	try:
		shutil.copy(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir1+'\\'+subdir2
	except IOError:
		print "Wrong path provided because the bed file does not exist."


# Get Genbank file
def genbank_file(accession_number):
	from subprocess import call
	cmd = ["python","Genomes_from_Genbank.py", accession_number]
	print 'calling: ' + ".".join(cmd)
	call(cmd)

# Run the genbank_to_bed12.py script
#import genbank_to_bed12
def accession_to_genome(genome):
	from subprocess import call
	cmd = ["python","genbank_to_bed12.py", accession_number]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	print 'accession_to_genome function running for ' + accession_number

# 5. Automate conversion from bed file to bigbed file on cheaha
# 	a. Create chrom.sizes file
def mk_chrom_sizes_file(genome,track_hub_directory,abrev):
	# 1. If known database, $ fetchChromSizes <db> > <db>.chrom.sizes
	# 2. If new database, determine the number of nucleotides in the chromosomes
	subdir = track_hub_directory
	filename = abrev+".chrom.sizes"
	dest_filepath = os.path.join(subdir, filename)
	file_contents = genome+"\t235646"
	chrom_sizes_file = open(dest_filepath,'w')
	chrom_sizes_file.write(file_contents)
	chrom_sizes_file.close
	'''
	# Save file to subdir	
	subdir = genome
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."
	'''
	''' 
	# To Fetch the Chrom Sizes When the Database Is Known: 
	from subprocess import call
	cmd = ["fetchChromSizes", accession_number]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	print 'accession_to_genome function running for ' + accession_number
	'''
# 	b. Sort the bed file: $ sort -k1,1 -k2,2n unsorted.bed > sorted.bed
def sort_bed_file(genome,abrev):
	old_filename = genome+'.bed'
	new_filename = genome+'sorted.bed'
	from subprocess import call
	cmd = ["sort","-k1,1","-k2,2n", old_filename,'-o', new_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	#print 'sort_bed_file function running for ' + accession_number
	# move sorted file to subdir
	# subdir = abrev
	# filename = new_filename
	# dest_filepath = os.path.join(subdir, filename)
	# try:
	# 	shutil.copyfile(filename,dest_filepath)
	# 	print filename+" moved to subdirectory "+ subdir
	# except IOError:
	# 	print "Wrong path provided because the bed file does not exist."

# 	c. Convert bed to bigbed
# 		1. $ module load ngs-ccts/ucsc_kent/2014-03-05
# 		2. $ bedToBigBed in.bed hg19.chrom.sizes out.bb
# 	d. Save bigbed file in directory
def bedToBigBed(genome):
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

# 6. Create the hub.txt file and save in directory on server
def mk_hub_txt_file(genome,track_hub_directory):
	hub_name = 'hub.txt'
	hub_short_label= 'Viruses-HH5(HCMV)'
	hub_long_label= 'Track hub for HCMV strain from RefSeq and GenBank'
	genomes_filelist= 'genomes.txt'
	email_address = 'hh5trackhub@vo.uabgrid.uab.edu'
	descriptionUrl = 'description.html'
	
	line1 = ' '.join(['hub',hub_name])
	line2 = ' '.join(['shortLabel', hub_short_label])
	line3 = ' '.join(['long_label', hub_long_label])
	line4 = ' '.join(['genomesFile', genomes_filelist])
	line5 = ' '.join(['email', email_address])
	line6 = ' '.join(['descriptionUrl', descriptionUrl])
	hub_list_of_lines = [line1,line2,line3,line4,line5,line6]
	filename = hub_name
	hub_file_open = open(filename,'w')
	hub_file_contents = '\n'.join([line1,line2,line3,line4,line5,line6])
	#print hub_file_contents
	for line in hub_list_of_lines:
		hub_file_open.writelines(line+'\n')
	hub_file_open.close()

	# Save file to subdir
	subdir = track_hub_directory
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

def mk_descriptionUrl_file(genome,track_hub_directory):
	filename = 'description.html'
	description_str = """<html>\n<body>\n<h1>hello world</h1>\nThis is HCMV (HH5) strain Merlin, for starters.\n\nSupport contact: hh5trackhub@vo.uabgrid.uab.edu\n</body>\n</html>"""
	#print description_str
	description_file = open(filename,'w')
	description_file.write(description_str)
	description_file.close

	# Save file to subdir
	subdir = track_hub_directory
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

def mk_genomes_file(genome,track_hub_directory):
	filename = 'genomes.txt'
	# genomes_str = 'genome '+abrev+'\ntrackDb '+abrev+'/trackDb.txt\ntwoBitPath '+abrev+'/NC_006273v2.2bit\ngroups '+abrev+'/groups.txt\ndescription '+genome+'\norganism HH5 strain Merlin\ndefaultPos NC_006273v2:1-235646\norderKey 100\nscientificName Human herpesvirus 5\nhtmlPath description.html'
	# genomes_file = open(filename,'w')
	# genomes_file.write(genomes_str)
	# genomes_file.close
	line1 = 'genome '+abrev
	line2 = 'trackDb '+abrev+'/trackDb.txt'
	line3 = 'twoBitPath '+abrev+'/NC_006273v2.2bit'
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

def mk_groups_file(genome,abrev):
	filename = 'groups.txt'
	groups_str = '''name gene\nlabel Genes and Gene Predictions\npriority 10\ndefaultIsClosed 1'''
	groups_file = open(filename,'w')
	groups_file.write(groups_str)
	groups_file.close

	# Save file to subdir
	# subdir = abrev
	# dest_filepath = os.path.join(subdir, filename)
	# try:
	# 	shutil.move(filename,dest_filepath)
	# 	print filename+" moved to subdirectory "+ subdir
	# except IOError:
	# 	print "Wrong path provided."

'''
8. Create the trackDb.txt files
	a. Minimum Requirements: 
		track track_name
		bigDataUrl track_data_URL
		shortLabel short_label
		longLabel long_label
		type track_type
	b. Track: symbolic name of track 
		must be unique
		must be the first entry in the trackDb.txt file
	c. BigDataUrl: the file name, path, or Web location of the track's data file
		Can be URL
		path relative to the trackDb.txt file is preferable because changes to the server name 
		won’t affect relative path, but will render the URL ineffective
	d. shortLabel: short name for track appears in left column of the genome browswer
		Max length 17
	e. longLabel
	“type bigBed #” where # is the number of features or columns in the bed file
9. Create track description files for each track and save in subdirectory
	Purpose: provide information about research to build credibility of data 
	and help people decide whether or not to use data
'''
def mktrackDb_file(genome,abrev):
	filename = "trackDb.txt"
	trackDb_str = "track "+abrev+"_refseq_mrna\nbigDataUrl "+genome+".bb\nshortLabel RefSeq Transcripts\nlongLabel RefSeq transcripts(mRNA)\ncolorByStrand 150,100,30 230,170,40\ncolor 150,100,30\naltColor 230,170,40\ntype bigBed 12\ngroup genes\nsearchIndex name"
	trackDb_file = open(filename,'w')
	trackDb_file.write(trackDb_str)
	trackDb_file.close

	# Save file to subdir
	# subdir = abrev
	# dest_filepath = os.path.join(subdir, filename)
	# try:
	# 	shutil.move(filename,dest_filepath)
	# 	print filename+" moved to subdirectory "+ subdir
	# except IOError:
	# 	print "Wrong path provided."

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
	# fasta file must be in same directory to call faToTwoBit below

	# Save file to subdir
	# subdir = abrev
	# dest_filepath = os.path.join(subdir, file_name)
	# try:
	# 	shutil.move(file_name,dest_filepath)
	# 	print file_name+" moved to subdirectory "+ subdir
	# except IOError:
	# 	print "Wrong path provided."

# add function to create 2bit file from fasta
def fasta_to_2bit(genome,abrev):
	from subprocess import call
	fasta_file = genome+'.fna'
	output_2bit_filename = genome+'.2bit'
	cmd = ["faToTwoBit",fasta_file, output_2bit_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)

def move_files(accession_number, abrev):
	gbk_file_name = accession_number+'.gbk'

	subdir = abrev
	dest_filepath = os.path.join(subdir, gbk_file_name)
	try:
		shutil.move(file_name,dest_filepath)
		print file_name+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

if __name__ == '__main__':
	mkdir_p(track_hub_directory)

	#genbank_file(accession_number)
	#copy_file(accession_number+'.gbk',track_hub_directory,abrev)
	#accession_to_genome(accession_number)
	#copy_file(genome+".bed",track_hub_directory,abrev)
	
	sort_bed_file(genome,abrev)
	copy_file(genome+"sorted.bed",track_hub_directory,abrev)
	mk_chrom_sizes_file(genome,track_hub_directory,abrev)
	bedToBigBed(genome)
	copy_file(genome+'.bb',track_hub_directory,abrev)
	#mk_hub_txt_file(genome,track_hub_directory)
	#mk_descriptionUrl_file(genome,track_hub_directory)
	mk_genomes_file(genome,track_hub_directory)
	copy_file(accession_number+'.gbk',track_hub_directory,abrev)
	change_dir(track_hub_directory)
	mkdir_p(abrev)
	change_dir(abrev)
	# check locations of file created in functions below
	mk_groups_file(genome,abrev)
	mktrackDb_file(genome,abrev)
	fasta_file(accession_number,abrev)
	fasta_to_2bit(genome,abrev)


# in project 0 folder on cheaha2
# 	$ cd project0
# load and activate biopython 
# 	$ module load Anaconda2/4.2.0
# 	$ source activate py27_biopython
# load Kent module
# 	$ module load Kent_tools/340

# run program with
# $  python accession_to_genome_w_exonerate.py NC_006273.2 Merlin_exon_v1
