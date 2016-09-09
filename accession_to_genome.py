#!/usr/local/bin/python
# coding: latin-1
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

# Use argv to obtain input of accession # from command line
from sys import argv 

accession_number = argv[1] # Upacks argv-> assigned to 1 variable you can work with


# Create directory, named by accession #, to place relevant files
from os import path
import errno    
import os
path_str = './'+str(accession_number)
def mkdir_p(path_str):
    try:
        os.makedirs(path_str)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path_str):
            pass
        else:
            raise

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
def accession_to_genome(accession_number):
	from subprocess import call
	cmd = ["python","genbank_to_bed12.py", "NC_006273.2"]
	print 'calling: ' + ".".join(cmd)
	call(cmd)
	print 'accession_to_genome function running for ' + accession_number


if __name__ == '__main__':
	mkdir_p(accession_number)
	#change_dir(accession_number)
	accession_to_genome(accession_number)
	

