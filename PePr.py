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

	peakfilename = opt.name+"__peak_shift-"
	for chip in readData.chip_filename_list:
		peakfilename += str(readData.shiftSize[chip])+'.'
	peakfilename = peakfilename[:-1]+"_window-"+str(windowsize)+"p-value-"+str(opt.threshold)
	sigTests.negative_binomial(readData,peakfilename,opt.peaktype,opt.remove_short_fragment,opt.narrow_peak_width,opt.threshold,windowsize)
#	area_disp.cal_shift_size_for_peak(readData,shiftSize)
	
if __name__ == '__main__':
	try: main(sys.argv)
	except KeyboardInterrupt: 
		print "user interrupted me"
		
		

