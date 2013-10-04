#!/usr/bin/env python26

import sys,re
import logging

root_logger = logging.getLogger("")
info = root_logger.info

def sam_parse(filename_list, chr_list):
	#parsing sam format files

	data_dict = {} #store the position and strand data for each file, seperated by chromosomes.
	infile = open(filename_list[0], 'r')

	for line in infile:
		words = line.strip().split()
		if not words[0].startswith("@"):
			length = len(words[9])
			break

	infile.close()
	info("the lenghth of the reads is %s", length)

	for chr in chr_list:
		data_dict[chr] = {}

		for file in filename_list:
			data_dict[chr][file] = {}

			for strand in ['f', 'r']:
				data_dict[chr][file][strand] = []

	for filename in filename_list:
		infile = open(filename, 'r')
		info("retrieving reads from file : %s", filename)
		for line in infile: 
			if not line.startswith("@"):
				break
		for line in infile:
			words = line.strip().split()
			words[1] = int(words[1])

			if not (words[1] & 0x0004):
				chr = words[2]
				pos = int(words[3])-1

				if not (words[1] & 0x16):
		#			print chr+':'+str(pos)+"-------"+'f'
					try:
						data_dict[chr][filename]['f'].append(pos)
					except KeyError: pass
				else:
		#			print chr+':'+str(pos)+"-------"+'r'
					try:
						data_dict[chr][filename]['r'].append(pos + length+1)
					except KeyError: pass
		infile.close()
	return data_dict,length

def bowtie_parse(filename_list, chr_list):
	#parsing bowtie format files
	data_dict = {} #store the position and strand data for each file, seperated by chromosomes.
	infile = open(filename_list[0], 'r')
	line = infile.readline()
	words = line.strip().split()
	length = len(words[4])
	info("the lenghth of the reads is %s", length)
	infile.close()

	for chr in chr_list:
		data_dict[chr] = {}

		for file in filename_list:
			data_dict[chr][file] = {}

			for strand in ['f', 'r']:
				data_dict[chr][file][strand] = []

	for filename in filename_list:
		infile = open(filename, 'r')
		info("retrieving reads from file : %s", filename)

		for line in infile:
			line = line.strip().split()
			chr = line[2]
			strand = line[1]
			pos = int(line[3])

			if strand == '+':
				try: data_dict[chr][filename]['f'].append(pos)
				except KeyError: pass
			elif strand == '-':
				try: data_dict[chr][filename]['r'].append(pos + length)
				except KeyError: pass
			else:
				print("strand error")
		infile.close()
	return data_dict,length

def bed_parse(filename_list,chr_list):
	#parsing BED format files
	data_dict={} #store the position and strand data for each file, seperated by chromosomes.
	for chr in chr_list:
		data_dict[chr]={}
		for file in filename_list:
			data_dict[chr][file]={}
			for strand in ['f','r']:
				data_dict[chr][file][strand]=[]

	for filename in filename_list:
		infile = open(filename,'r')
		info("retrieving reads from file: %s",filename)
		for line in infile:
			line = line.strip().split()
			chr=line[0]

			strand = line[5]
			if strand == '+':
				pos = int(line[1])
				try: data_dict[chr][filename]['f'].append(pos)
				except KeyError: pass
			elif strand == '-':
				pos = int(line[2])-1
				try: data_dict[chr][filename]['r'].append(pos)
				except KeyError: pass

			else:
				print("strand error")
		infile.close()
	return data_dict,int(line[2])-int(line[1]) #readLength

def process_illumina_match(align,mismatch,length,format):
	words = align.split(',')
	if format == "multi":
		for match in words:
			if match.startswith('c'):
				str,match =match.split(':')
				chr = re.match(r'(chr(\d\d?|\w))',str)
				chr = chr.group(0)
			if match.endswith(mismatch):
				strand = match[-2]
				pos = int(match[0:-2])
				if strand =='R':
					pos = pos+length
				return chr,strand,pos
	elif format == "extended":
		for match in words:
			if match.startswith('c'):
				chr_part, match = match.split(':')
				chr = re.match(r'(chr(\d\d?|\w))', chr_part)
				chr = chr.group(0)
				m1 = re.search(r'([F|R])', match)
				strand = m1.group(1)
				pos = int(re.sub(r'[F|R].*', '', match))

				if strand == 'R':
					pos = pos + length

				match = re.sub(r'(\^\w+\$)', '', match)
				m2 = re.findall(r'[A-Z]', match)

				if len(m2) - 1 == int(mismatch):
					return chr, strand, pos


def eland_parse(filename_list,chr_list,format="default"):
	#parsing eland format files
	data_dict={} #store the position and strand data for each file, seperated by chromosomes.
	for chr in chr_list:
		data_dict[chr]={}
		for file in filename_list:
			data_dict[chr][file]={}
			for strand in ['f','r']:
				data_dict[chr][file][strand]=[]
## initialize the data structure
	infile = open(filename_list[0],'r')
	line = infile.readline()
	words=line.split()
	length=len(words[1])
	info ("the length of the reads is %s",length)
	infile.close()

	for filename in filename_list:
		infile = open(filename,'r')
		info("retrieving reads from file:%s",filename)
		for line in infile:
			words = line.strip().split()
			mismatch = -1
			
			if len(words) == 4 and ':' in words[2]:
				num_mismatch = words[2].split(':')
				
				if num_mismatch[0]== '1':
					mismatch = '0'
				elif (num_mismatch[0]=='0') & (num_mismatch[1] == '1'):
					mismatch = '1'
				elif (num_mismatch[0]=='0') & (num_mismatch[1]=='0') & (num_mismatch[2] == '1'):
					mismatch = '2'
				if mismatch !=-1:
					try:
						chr,strand,pos = process_illumina_match(words[3],mismatch,length,format)
						if strand == 'F':
							try: data_dict[chr][filename]['f'].append(pos)
							except KeyError:pass
						elif strand == 'R':
							try: data_dict[chr][filename]['r'].append(pos)
							except KeyError:pass
					except TypeError:pass
		infile.close()
	return data_dict,length



def parse(readData,fileFormat):
	''' read in the supported file format and make it into the data file that we need'''
	info("begin parsing files...")
	data_dict = {}
	if fileFormat == "bed":
		data_dict,readLength = bed_parse(readData.filename_list,readData.chr_list)
	elif fileFormat == "eland_multi":
		data_dict,readLength = eland_parse(readData.filename_list,readData.chr_list,format="multi")
	elif fileFormat == "eland_extended":
		data_dict,readLength = eland_parse(readData.filename_list,readData.chr_list,format="extended")
	elif fileFormat == "bowtie":
		data_dict,readLength = bowtie_parse(readData.filename_list,readData.chr_list)
	elif fileFormat == "sam":
		data_dict,readLength = sam_parse(readData.filename_list,readData.chr_list)
	else: 
		pass

	#### here can add other file format support ####
	
	#fill in the data structures. 
	readData.data_dict_by_strands = data_dict
	readData.readLength = readLength
	readData.cal_read_total()
	readData.refine_chr_list()
	readData.cal_genomeSize()
	info( "leaving file_parser.py")

def main():

	filename_list = sys.argv[1:]
	data_dict = {}
	chr_list = []
	data_dict,chr_list = parse("bed", filename_list)
	



if __name__ == '__main__':
	main()
