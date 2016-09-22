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

# Get Genbank file
def genbank_file(accession_number):
	'''
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
	'''
	from subprocess import call
	cmd = ["python","Genomes_from_Genbank.py", accession_number]
	print 'calling: ' + ".".join(cmd)
	call(cmd)
	#print 'genbank_file function running for ' + accession_number


# if accession number has a period, indicating version number, replace with v.
if '.' in accession_number:
	path_str='v'.join(accession_number.split('.'))
else:
	path_str = str(accession_number)
#print path_str

def mkdir_p(path_str):
	try:
		os.makedirs(path_str)
	except WindowsError as exc:# Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path_str):
			pass
		else:
			raise
	print "New directory, "+path_str+", created in the current directory."
	subdir = path_str
	filename = path_str+".bed"
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.copyfile(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided because the bed file does not exist."

def change_dir(path_str):
	try:
		os.chdir(path_str)
	except WindowsError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path_str):
			pass
		else:
			raise

# Run the genbank_to_bed12.py script
#import genbank_to_bed12
def accession_to_genome(path_str):
	from subprocess import call
	cmd = ["python","genbank_to_bed12.py", accession_number]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	print 'accession_to_genome function running for ' + accession_number

# 5. Automate conversion from bed file to bigbed file on cheaha
# 	a. Create chrom.sizes file
def mk_chrom_sizes_file():
	# 1. If known database, $ fetchChromSizes <db> > <db>.chrom.sizes
	# 2. If new database, determine the number of nucleotides in the chromosomes
	filename = "hh5Merlin2.chrom.sizes"
	file_contents = "NC_006273v2	235646"
	chrom_sizes_file = open(filename,'w')
	chrom_sizes_file.write(file_contents)
	chrom_sizes_file.close

	# Save file to subdir
	subdir = path_str
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.copyfile(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."
	'''
	cmd = ["fetchChromSizes", accession_number]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	print 'accession_to_genome function running for ' + accession_number
	'''
# 	b. Sort the bed file: $ sort -k1,1 -k2,2n unsorted.bed > sorted.bed
def sort_bed_file(path_str):
	old_filename = path_str+'.bed'
	new_filename = path_str+'sorted.bed'
	from subprocess import call
	cmd = ["sort","-k1,1","-k2,2n", old_filename,'-o', new_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	#print 'sort_bed_file function running for ' + accession_number
	# move sorted file to subdir
	subdir = path_str
	filename = new_filename
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.copyfile(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided because the bed file does not exist."

# 	c. Convert bed to bigbed
# 		1. $ module load ngs-ccts/ucsc_kent/2014-03-05
# 		2. $ bedToBigBed in.bed hg19.chrom.sizes out.bb
# 	d. Save bigbed file in directory
def bedToBigBed(path_str):
	sorted_bed_file = path_str+'sorted.bed'
	output_filename = path_str+'.bb'
	from subprocess import call
	cmd = ["bedToBigBed",sorted_bed_file,"hh5Merlin2.chrom.sizes", output_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)

# 6. Create the hub.txt file and save in directory on server
def mk_hub_txt_file(path_str):
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
	subdir = path_str
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

def mk_descriptionUrl_file(path_str):
	filename = 'description.html'
	description_str = """<html>\n<body>\n<h1>hello world</h1>\nThis is HCMV (HH5) strain Merlin and strain BE/7/2011, for starters.\n\nSupport contact: hh5trackhub@vo.uabgrid.uab.edu\n</body>\n</html>"""
	#print description_str
	description_file = open(filename,'w')
	description_file.write(description_str)
	description_file.close

	# Save file to subdir
	subdir = path_str
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

def mk_genomes_file(path_str):
	filename = 'genomes.txt'
	genomes_str = '''genome hh5Merlin2\ntrackDb hh5Merlin2/trackDb.txt\ntwoBitPath hh5Merlin2/NC_006273v2.2bit\ngroups hh5Merlin2/groups.txt\ndescription NC_006273v2\norganism HH5 strain Merlin\ndefaultPos NC_006273v2:1-235646\norderKey 100\nscientificName Human herpesvirus 5\nhtmlPath description.html'''
	#print description_str
	genomes_file = open(filename,'w')
	genomes_file.write(genomes_str)
	genomes_file.close

	# Save file to subdir
	subdir = path_str
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

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
def mktrackDb_file(path_str):
	filename = "trackDb.txt"
	trackDb_str = "track hh5Merlin2_refseq_mrna\nbigDataUrl NC_006273v2.mrna.bb\nshortLabel RefSeq Transcripts\nlongLabel RefSeq transcripts(mRNA)\ncolorByStrand 150,100,30 230,170,40\ncolor 150,100,30\n\naltColor 230,170,40\ntype bigBed 12\ngroup genes\n\ntrack hh5Merlin2_refseq_loci\nbigDataUrl NC_006273v2_refseq_loci.bb\nshortLabel RefSeq loci\nlongLabel RefSeq loci\ncolorByStrand 150,100,30 230,170,40\ncolor 150,100,30\naltColor 230,170,40\ntype bigBed 6\nsearchIndex name\ngroup genes"
	trackDb_file = open(filename,'w')
	trackDb_file.write(trackDb_str)
	trackDb_file.close

	# Save file to subdir
	subdir = path_str
	dest_filepath = os.path.join(subdir, filename)
	try:
		shutil.move(filename,dest_filepath)
		print filename+" moved to subdirectory "+ subdir
	except IOError:
		print "Wrong path provided."

if __name__ == '__main__':
	genbank_file(accession_number)
	accession_to_genome(accession_number)
	mkdir_p(path_str)
	sort_bed_file(path_str)
	mk_chrom_sizes_file()
	bedToBigBed(path_str)
	mk_hub_txt_file(path_str)
	mk_descriptionUrl_file(path_str)
	mk_genomes_file(path_str)
	mktrackDb_file(path_str)
# $ python accession_to_genome.py NC_006273.2
