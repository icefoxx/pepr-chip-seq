import sys
import misc
import logging
import logConfig
import optParser
import fileParser
import numpy


root_logger =  logging.getLogger("")
info = root_logger.info
debug = root_logger.debug
warning = root_logger.warning
shift_logger = logging.getLogger("shiftSizeEst")

def shift_size_per_chrom(forward, reverse,file=-1):
	#estimate the shift size for each chromosome separately. 

	shift_max = 200
	
	max_overlapping =  len(list(set(forward) & set(reverse)))
	optimum_shift = 0
	
	# iterate from 0 offset to shift_max, to find which offset produce the maximal overlap.  
	for offset in xrange(1,shift_max):
		
		overlapping = set([x+offset for x in forward]) & set([y-offset for y in reverse])
		overlapping = len(list(overlapping))

		if overlapping > max_overlapping:
			max_overlapping = overlapping
			optimum_shift = offset
		if file !=-1:
			file.write( str(offset)+'\t'+str(overlapping) + '\n')
	return optimum_shift

def estimate_shift_size(readData,opt):
	if opt.shiftSize != "-1": #if shift size provided by the user, then skip it.
		shift_list = opt.shiftSize.split(',')
		if len(shift_list)== 1:
			for name in readData.filename_list: 
				readData.shiftSize[name] = int(shift_list[0])
				info("%-10s %s",name,shift_list[0])
		else: 
			if opt.diff: 
				for idx,name in enumerate(readData.filename_list):
					readData.shiftSize[name] = int(shift_list[idx])
					info("%-10s %s",name,shift_list[idx])
			else: 
				for idx,chip in enumerate(readData.chip_filename_list):
					readData.shiftSize[chip] = int(shift_list[idx])
					info("%-10s %s",chip,shift_list[idx])
				for control in readData.control_filename_list:
					readData.shiftSize[control] = sum(shift_list)/len(shift_list)
					info("%-10s %s",control,sum(shift_list)/len(shift_list))
		return 
	info("begin estimating the shift size...")
	if not opt.diff: #if it is test samples versus the control 
		for chip_filename in readData.chip_filename_list:
			info("estimating for %-10s", chip_filename)
			shift_list = []

			for count,chr in enumerate(readData.chr_list):
				if count ==5:   #estimate the shift size for five chromosomes and take the median of these as the estimator. 
					break
				forward = readData.data_dict_by_strands[chr][chip_filename]['f']
				reverse = readData.data_dict_by_strands[chr][chip_filename]['r']
				shift = shift_size_per_chrom(forward,reverse)
				info("%-10s %d", chr, shift)
				shift_list.append(shift)
			shift_median = misc.median(shift_list[:])
			info("%-10s %d",chip_filename,shift_median)
			readData.shiftSize[chip_filename] = shift_median
		con_shift_size = sum(readData.shiftSize.values())/len(readData.shiftSize.values())
		for control_filename in readData.control_filename_list: 
			readData.shiftSize[control_filename] = con_shift_size
			
	else: #else, if we're calling differential binding. 
		for filename in readData.filename_list: 
			info("estimating for %-10s", filename)
			shift_list = []
			for count,chr in enumerate(readData.chr_list):
				if count==5: 
					break
				forward = readData.data_dict_by_strands[chr][filename]['f']
				reverse = readData.data_dict_by_strands[chr][filename]['r']
				shift = shift_size_per_chrom(forward,reverse)
				info("%-10s %d", chr, shift)
				shift_list.append(shift)
			shift_median = misc.median(shift_list[:])
			info("%-10s %d",filename,shift_median)
			readData.shiftSize[filename] = shift_median
	return 

	
def shift_reads(readData): #shift the reads according the specific shift size estimated or provided by user
	readData.data_dict = {}
	for chr in readData.chr_list:
		debug("shifting the reads for %s...",chr)
		readData.data_dict[chr] ={}
		for file in readData.filename_list:
			forward=readData.data_dict_by_strands[chr][file]['f']
			reverse=readData.data_dict_by_strands[chr][file]['r']
			readData.data_dict_by_strands[chr][file]['f'] = numpy.array(readData.data_dict_by_strands[chr][file]['f'])
			readData.data_dict_by_strands[chr][file]['r'] = numpy.array(readData.data_dict_by_strands[chr][file]['r'])
			readData.data_dict[chr][file]=[x+readData.shiftSize[file] for x in forward]+[y-readData.shiftSize[file] for y in reverse]

#	return data_dict			



def main(argv):
## for the purpose of estimating the shift size only

	# initialize the logger
	root_logger = logging.getLogger("")
	debug = root_logger.debug
	info = root_logger.info
	## performing the option parser
	opt = optParser.optParser(argv)
	optParser.validateOpt(opt)
	readData = opt.read_data
	## read in the data
	fileParser.parse(readData, opt.fileFormat)

	## remove the redundant reads
	if (opt.remove_redundant):
		readData.remove_redundant_reads()

	## shiftSize estimation and shifting reads
	shiftSize = estimate_shift_size(readData,opt)
	info (" The shiftSize is %s", shiftSize)

if __name__ == '__main__':
	main(sys.argv)
