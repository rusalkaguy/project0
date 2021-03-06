from __future__ import division
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
import pprint
pp = pprint.PrettyPrinter(indent=0)



bed_file_name = argv[1] # Upacks argv-> assigned to 1 variable you can work with
abrev = argv[2] # abrev is the the abreviated name of the chrom.sizes file and subdirectory, which should be SR10-01
genome = abrev
track_hub_directory = 'track_hub'
preferred_genome_accession = 'NC_'
# Read in bed file
bed_file=open(bed_file_name,'r')
line_list =bed_file.readlines()
col_list = line_list[1].split('\t')
# accession_number = col_list[0]
# genome = accession_number
#print line_list
num_lines = len(line_list)
#print num_lines
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
bed_file.close()
for line in output:
	processed_bed_file.write(line)
processed_bed_file.close()



def mkdir_p(abrev):
	try:
		os.makedirs(abrev)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(abrev):
			pass
		else:
			raise
	print "New directory, "+abrev+", created in the current directory."

def sort_bed_file(processed_bed_file_name, genome):
	old_filename = processed_bed_file_name
	new_filename = genome+'sorted.bed'
	from subprocess import call
	cmd = ['sort','-k1,1','-k2,2n', old_filename,'-o', new_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	
def create_normalized_bed_files(sorted_bed_filename,genome):
	# Normalize score from 0-1000 per gene per chromosome
	sorted_bed_filename = genome+'sorted.bed'
	normalized_bed_filename = genome+'.bed'
	best_normalized_bed_filename = genome+'best.bed'
	# Read in sorted bed file
	sorted_bed_file = open(sorted_bed_filename,'r')
	line_list = sorted_bed_file.readlines()
	contigs = {}
	genes = {}
	all_bed_output = []
	linenum = 0
	for line in line_list:
		linenum = int(linenum+1) #one-based numbering of file lines
		col_list = line.split('\t')
		#print col_list
		#print 'length of col_list:', len(col_list)
		#print line
		chrom = col_list[0]
		chromStart = col_list[1]
		chromEnd = col_list[2]
		#gene_name = col_list[3]
		gene_name_and_genome = col_list[3].split('.')
		gene_name = gene_name_and_genome[0]
		genome_accession = gene_name_and_genome[1]
		#print gene_name +','+ genome_accession
		score = col_list[4]
		strand = col_list[5]
		thickStart = col_list[6]
		thickEnd = col_list[7]
		itemRgb = col_list[8]
		blockCount = col_list[9]
		blockSizes = col_list[10]
		blockStarts = col_list[11].replace("\n", "")

		# Make Dictionary of Contigs per chromosome
		# contigs = {chrom:{gene_name:score}}
		# genes = {gene_name:score}
		if chrom not in contigs.keys(): 
			# initialize genes dict because different gene_names per chrom
			genes = {}
			#genes = {gene_name:score}
			genes[gene_name] = { "max_score": int(score), "linenum": int(linenum)}
			contigs[chrom] = genes
		if gene_name not in genes.keys():
			#genes[gene_name] = score
			#contigs[chrom] = genes
			contigs[chrom][gene_name] = { "max_score": int(score), "linenum": int(linenum)}		
		else:
			#genes = {gene_name,max(score,contigs[chrom][gene_name])}
			#contigs[chrom] = genes
			#contigs[chrom][gene_name] = max(int(score),int(contigs[chrom][gene_name]))
			if int(score)>int(contigs[chrom][gene_name]["max_score"]):
				contigs[chrom][gene_name]["max_score"] = int(score)
				contigs[chrom][gene_name]["linenum"] = int(linenum)
			if int(score)==int(contigs[chrom][gene_name]["max_score"]) and preferred_genome_accession == genome_accession[0:2]:
				contigs[chrom][gene_name]["linenum"] = int(linenum)

		orig_bed_line = '\t'.join([chrom,chromStart,chromEnd,gene_name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts,genome_accession])
		all_bed_output.append(orig_bed_line+'\n')
	#pp.pprint(contigs)
	#print all_bed_output
	sorted_bed_file.close()
	# write all_bed_output to bed file with 12+1 cols
	processed_bed_file = open(normalized_bed_filename,'w')
	for line in all_bed_output:
		processed_bed_file.write(line)
	processed_bed_file.close()
	# Open newly written file to normalize scores
	processed_bed_file = open(normalized_bed_filename,'r')
	processed_line_list = processed_bed_file.readlines()
	all_bed_output = []
	best_bed_output = []
	linenum = 0
	#print processed_line_list
	for line in processed_line_list:
		linenum = linenum+1 #one-based numbering of file lines
		col_list = line.split('\t')
		chrom = col_list[0]
		chromStart = col_list[1]
		chromEnd = col_list[2]
		gene_name = col_list[3]
		score = col_list[4]
		strand = col_list[5]
		thickStart = col_list[6]
		thickEnd = col_list[7]
		itemRgb = col_list[8]
		blockCount = col_list[9]
		blockSizes = col_list[10]
		blockStarts = col_list[11]
		genome_accession = col_list[12]
		# Normalize and scale the score
		normalized_score = int(round((int(score)/int(contigs[chrom][gene_name]["max_score"]))*1000))
		scaled_norm_score = str(int(round(max(20*(normalized_score-950),0))))
		# Reconstructed bed line with scaled, normalized score
		normalized_bed_line = '\t'.join([chrom,chromStart,chromEnd,gene_name,scaled_norm_score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts,genome_accession])
		all_bed_output.append(normalized_bed_line)
		
		# Write bed output with only the best matches (scores)

		if int(linenum) ==int(contigs[chrom][gene_name]["linenum"]):
			best_bed_output.append(normalized_bed_line)
			
	#pp.pprint(best_bed_output)
	processed_bed_file.close()
	#print output
	# Open file, write all scaled and normalized bed lines, and close file
	all_normalized_bed_file = open(normalized_bed_filename,'w')
	for line in all_bed_output:
		all_normalized_bed_file.write(line)
	all_normalized_bed_file.close()
	# Open file, write only best scaled and normalized matches, and close file
	best_normalized_bed_file = open(best_normalized_bed_filename, 'w')
	for line in best_bed_output:
		best_normalized_bed_file.write(line)
	best_normalized_bed_file.close()

	bed_filenames_list= [normalized_bed_filename,best_normalized_bed_filename]
	return bed_filenames_list

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
	# 3. If multiple chromosomes or contigs, count the number of nucleotides in each

	subdir = track_hub_directory
	filename = abrev+".chrom.sizes"
	fasta_filename = abrev+'.fa'
	dest_filepath = os.path.join(subdir, fasta_filename)
	# for 1 chromosome:
	#file_contents = genome+"\t235646"
	# for more than 1 chromosome/contig:
	# count the number of nucleotides in each chrom/contig from fasta file
	import sys
	from Bio import SeqIO

	FastaFile = open(fasta_filename, 'rU')
	output = []
	for rec in SeqIO.parse(FastaFile, 'fasta'):
		name = rec.id
		seq = rec.seq
		seqLen = len(rec)
		#print name, seqLen
		chrom_sizes_line = '\t'.join([name,str(seqLen)])
		output.append(chrom_sizes_line+'\n')
	#print output
	FastaFile.close()
	# write output to file in subdir
	chrom_sizes_file = open(dest_filepath,'w')
	for line in output:
		chrom_sizes_file.write(line)
	chrom_sizes_file.close()
	#chrom_sizes_file.write(output)
	print 'wrote '+filename+' to '+dest_filepath
	# write output to file in current directory
	chrom_sizes_file = open(filename,'w')
	for line in output:
		chrom_sizes_file.write(line)
	chrom_sizes_file.close()
	print 'wrote '+filename+' to current directory'


def bedToBigBed(input_bed_filename,abrev):
	from subprocess import call
	#normalized_bed_file = genome+'normalized.bed'
	#output_filename = genome+'.bb'
	output_filename = input_bed_filename.split('.')[0]+'.bb'
	chrom_sizes_file = abrev+'.chrom.sizes'
	as_file = '-as=gene_prediction.fields.as'
	# bedToBigBed -extraIndex=name -type=bed12 fileinsorted.bed chrom.sizes outputfile.bb
	#index = '-extraIndex=name'
	#bed_type = '-type=bed12'
	bed_type = '-type=bed12+1'
	index = '-extraIndex=name,genomeAccession'
	cmd = ["bedToBigBed",as_file,index,bed_type, input_bed_filename, chrom_sizes_file, output_filename]
	#cmd = ["bedToBigBed",index,bed_type, normalized_bed_file, chrom_sizes_file, output_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)
	return output_filename
# change entries for this genome and track hub
def mk_hub_txt_file(genome,track_hub_directory):
	hub_name = 'hub.txt'
	hub_short_label= 'Viruses-'
	hub_long_label= 'Track hub for contigs of all genomes from Vicuna'
	genomes_filename= genome+'.genomes.txt'
	email_address = 'bheater.uab.edu'
	descriptionUrl = 'description.html'
	
	line1 = ' '.join(['hub',hub_name])
	line2 = ' '.join(['shortLabel', hub_short_label])
	line3 = ' '.join(['longLabel', hub_long_label])
	line4 = ' '.join(['genomesFile', genomes_filename])
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
	description_str = """<html>\n<body>\n<h1>hello world</h1>\nEnter information about this strain of virus here.\n\nSupport contact: bheater.uab.edu\n</body>\n</html>"""
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
	# read in length of chrom 1 from chrom.sizes file SR10-01.chrom.sizes
	chrom_sizes_file = abrev+'.chrom.sizes'
	chrom_sizes_file_open = open(chrom_sizes_file,'r')
	first_line = chrom_sizes_file_open.readline()
	#
	chrom_sizes_file_open.close()
	filename = genome+'.genomes.txt'
	line1 = 'genome '+abrev
	line2 = 'trackDb '+abrev+'/trackDb.txt'
	line3 = 'twoBitPath '+abrev+'/'+genome+'.2bit'
	line4 = 'groups '+abrev+'/groups.txt'
	line5 = 'description '+genome
	line6 = 'organism Virus'
	line7 = 'defaultPos '+genome+':1-235646'
	line8 = 'orderKey 100'
	#line9 = 'scientificName Human herpesvirus 5'
	line10 = 'htmlPath description.html'
	genomes_list_of_lines = [line1,line2,line3,line4,line5,line6,line7,line8,line10]
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
	groups_file.close()

def mktrackDb_file(input_bigbed_filenames_list,abrev):
	trackDb_list = []
	track_filename = "trackDb.txt"
	visibility = "visibility 3"
	shortLabel = "shortLabel Virus Genes"
	longLabel = "longLabel Virus Genes with Best Scores"
	#colorByStrand= "colorByStrand 150,100,30 230,170,40\ncolor 150,100,30\naltColor 230,170,40\n
	#color = "color 0,60,120 " use different colors for genes, mrna, transcripts later on
	useScore="useScore 1"
	trackType = "type bigBed"
	group = "group genes"
	searchIndex = "searchIndex name,genomeAccession"
	for filename in input_bigbed_filenames_list:
		#print filename.split('.')[0] ==abrev +'best'
		#print abrev + 'best'
		#print filename.split('.')[0]
		if filename.split('.')[0] ==abrev +'best':
			track = "track "+'best_'+abrev+"_genes"
			longLabel = "longLabel Virus Genes with Best Scores"
		else: 
			track = "track "+'all_'+abrev+"_genes"
			longLabel = "longLabel All Virus Genes"
		bigDataUrl = "bigDataUrl "+filename
		trackDb_str = '\n'.join([track,visibility,bigDataUrl,shortLabel,longLabel,useScore,trackType,group,searchIndex])
		trackDb_list.append(trackDb_str+'\n')
	
	#pp.pprint(trackDb_list)
	#print "filename: " + track_filename
	trackDb_file = open(track_filename,'w')
	for line in trackDb_list:
		trackDb_file.write(line+'\n')
	trackDb_file.close()

def fasta_to_2bit(abrev,genome):
	from subprocess import call
	subdir = abrev
	output_2bit_filename = genome+'.2bit'
	#dest_filepath = os.path.join(subdir, output_2bit_filename)
	fasta_file = genome+'.fa'
	cmd = ["faToTwoBit",fasta_file, output_2bit_filename]
	print 'calling: ' + " ".join(cmd)
	call(cmd)

bigbed_filenames_list = []
if __name__ == '__main__':
	mkdir_p(track_hub_directory)
	sort_bed_file(processed_bed_file_name,genome)
	copy_file(genome+"sorted.bed",track_hub_directory,abrev)
	bed_filenames_list = create_normalized_bed_files(genome+'sorted.bed',genome)
	for bed_filename in bed_filenames_list:
		copy_file(bed_filename,track_hub_directory,abrev)
	mk_chrom_sizes_file(genome,track_hub_directory,abrev)
	#copy chrom_sizes file to track_hub directory
	try:
		shutil.copy(genome+".chrom.sizes",os.path.join(track_hub_directory, genome+".chrom.sizes"))
		print genome+".chrom.sizes moved to subdirectory "+ track_hub_directory
	except IOError:
		print "Wrong path provided because the bed file does not exist."

	for input_bed_filename in bed_filenames_list:
		ouput_bigbed_filename = bedToBigBed(input_bed_filename,abrev)
		bigbed_filenames_list.append(ouput_bigbed_filename)
	#print bigbed_filenames_list
	for bigbed_filename in bigbed_filenames_list:
		copy_file(bigbed_filename,track_hub_directory,abrev)

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
	mk_groups_file(genome,abrev)
	mktrackDb_file(bigbed_filenames_list,abrev)
	print '\n'

# in project 0 folder on cheaha2
# 	$ cd project0
# load and activate biopython 
# 	$ module load Anaconda2/4.2.0
# 	$ source activate py27_biopython
# load Kent module
# 	$ module load -i Kent_tools/340


# to run
# python bed_to_trackhub.py all_geno.SR10-01.vicuna.contig.score100.best1.bed SR10-01