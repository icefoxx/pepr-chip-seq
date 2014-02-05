#!/usr/bin/python26
import re, os, sys,time
import logging
import logConfig
import optParser
import fileParser
import shiftSize
import windowSize
import sigTests
import genomeData
import misc

def main(argv):
	# initialize the logger 
	root_logger = logging.getLogger("")
	debug = root_logger.debug
	info = root_logger.info
	## performing the option parser
	opt = optParser.optParser(argv)
	optParser.validateOpt(opt)
	readData = opt.read_data
	##1. read and parse the data
	fileParser.parse(readData, opt.fileFormat)

	##2. remove the redundant reads
	if (opt.remove_redundant):
		readData.remove_redundant_reads()	

	##3. shiftSize estimation and shifting reads
	shiftSize.estimate_shift_size(readData,opt)
	shiftSize.shift_reads(readData)
	##4. windowSize estimation and split reads into windows
	windowsize = windowSize.estimate_window_size(readData,opt)
	info (" The windowSize is %s", windowsize)
	windowSize.separate_exact_by_window(readData, windowsize) 

	# del readData.data_dict ##could reduce the memory usage
	if not opt.diff:
		swap = False
		peakfilename = opt.name+"__PePr_peaks.bed"
		sigTests.negative_binomial(readData,peakfilename,opt.peaktype,swap,opt.remove_short_fragment,opt.narrow_peak_width,opt.threshold,windowsize)
	else: 
		up_peakfilename = opt.name+"__PePr_up_peaks.bed"
		down_peakfilename = opt.name+"__PePr_down_peaks.bed"
		swap = False
		sigTests.negative_binomial(readData,up_peakfilename,opt.peaktype,swap,opt.remove_short_fragment,opt.narrow_peak_width,opt.threshold,windowsize)
		swap = True
		sigTests.negative_binomial(readData,down_peakfilename,opt.peaktype,swap,opt.remove_short_fragment,opt.narrow_peak_width,opt.threshold,windowsize)
	#write to a file that record the command and parameters. 	

	
if __name__ == '__main__':
	try: main(sys.argv)
	except KeyboardInterrupt: 
		print "user interrupted me"
		
		

