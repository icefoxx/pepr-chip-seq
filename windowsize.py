#!/usr/bin/env python

import sys,re,math
import time
import random
import numpy
import numpy.lib.stride_tricks as std
import misc
import logging

root_logger = logging.getLogger("")
debug = root_logger.debug
info = root_logger.info
error = root_logger.error

window_logger = logging.getLogger("windowSizeEst")

def get_window_size2(array,bin=20,iter=100):
	# this is the function that estimate the window size. 
        chr_len = array.size
        rowNum = chr_len/bin
        array_window = std.as_strided(array,(rowNum,bin),(bin*array.itemsize,1*array.itemsize))
        array_window = numpy.sum(array_window,1)

        peak_len_list = []
        for x in range(iter):
		peak = numpy.max(array_window)
		peak_idx = numpy.where(array_window == peak)[0][0]

                array_window[peak_idx]= -1 # set the peak to -1 so that it won't confuse the next iteration.

		
		if peak_idx == 0:  #check if the window reaches the left boundary 
			left_boundary_not_reached = False
		else:
			left_boundary_not_reached = True	
		if peak_idx == rowNum-1: #check if the window reaches of the right boundary 
			right_boundary_not_reached = False
		else: 
			right_boundary_not_reached = True
		i_l = 1
                i_r = 1 # be careful if anything equals to -1. though it's unlikely that it will happen.

                while left_boundary_not_reached:

                        while (array_window[peak_idx-i_l] >= 0.2*peak): # probably subtract background in future
                                array_window[peak_idx-i_l] = -1
				if peak_idx-i_l==0:
					left_boundary_not_reached = False
					break
                                i_l += 1
			if left_boundary_not_reached:
				if peak_idx-i_l==0:
					array_window[peak_idx-i_l] = -1
					left_boundary_not_reached = False
			
                        if left_boundary_not_reached: 
				if (array_window[peak_idx-i_l-1] >= 0.2*peak): #will continue counting if the next window(one window gap) has above the 20% mode reads. 
					array_window[peak_idx-i_l] = -1
                                	i_l += 1
                        	else:
                        	        break

                while right_boundary_not_reached:
                        while (array_window[peak_idx+i_r] >= 0.2*peak):
                                array_window[peak_idx+i_r] = -1
				if peak_idx+i_r ==rowNum-1:
					right_boundary_not_reached = False
					break
                                i_r += 1
			if right_boundary_not_reached:
				if peak_idx+i_r==rowNum-1:
					array_window[peak_idx+i_r] = -1
					right_boundary_not_reached = False

                        if right_boundary_not_reached:
				if (array_window[peak_idx+i_r+1] >=0.2*peak): #will continue counting if the next window(one window gap) has above the 20% mode reads.
					array_window[peak_idx+i_r] = -1
                                	i_r += 1
                        	else:
                                	break
                peak_len_list.append((1+i_l+i_r)*bin)

        return misc.median(peak_len_list)

def estimate_window_size(readData,opt):
	if (opt.windowSize != -1):
		# window size should be even number 
		if opt.windowSize % 2 ==1:
			info("warning: window size should be even number; adding 1 to it.")
			opt.windowSize = opt.windowSize + 1
		return opt.windowSize

	info( "begin window_size calculation...")
	window_size_list = []
	for chr in readData.chr_list:
		reads = []
		for chip_filename in readData.chip_filename_list:
			reads.extend(readData.data_dict[chr][chip_filename])
		coord_array = numpy.zeros(readData.chr_length_dict[chr],dtype=numpy.int16)
		for x in reads:
			try: coord_array[x] += 1
			except IndexError: pass #debug("coordinates out of range")
#		debug("the length of the %-5s array is %s", chr, coord_array.size)

		window = get_window_size2(coord_array)
		if window >1000: # limit the window size to be less than 1000 bases...
			window = 1000
		window_size_list.append(window)
		
		window_logger.debug("%-10s %s", chr, window)
	
	debug("length of windowsize list is "+str(len(window_size_list)))
	window_median = misc.median(window_size_list)
	# window size cannot be odd number since we need to use overlapping half windows. 
	if window_median % 2 == 1:
		window_median = window_median+1
	window_logger.debug("windowSize: %+10s",window_median)
	info ( "finishing windowSize calculation...")
	return window_median 


def separate_exact_by_window(readData,windowSize,normalize = "Large"):
	data_by_window_dict = {}
	index_dict = {}
	scale_dict = {}
	if normalize =="Large":	#scale up to the largest sample
		read_max = float(max(readData.read_total_per_file.values()))
		for filename in readData.filename_list:
			scale_dict[filename] = read_max / readData.read_total_per_file[filename]
			debug("The scaling up factor for %s is %s",filename,scale_dict[filename])
	elif normalize =="Small":#scale down to the smallest sample
		read_min = float(min(readData.read_total_per_file.values()))
		for filename in readData.filename_list:
			scale_dict[filename] = read_min / readData.read_total_per_file[filename]
			debug("The scaling up factor for %s is %s",filename,scale_dict[filename])
			

	move_size = windowSize/2

	for chr in readData.chr_list:
		info( "partitioning by window on \t%s...", chr)
		chr_len =readData.chr_length_dict[chr]
		rowNum = chr_len/move_size - 1

		for idx,filename in enumerate(readData.filename_list):
			array = numpy.zeros(rowNum,dtype=numpy.int32)
			for x in readData.data_dict[chr][filename]:
				try: array[x/move_size]+=1
				except IndexError:    pass  # debug("coordinates out of range"+str(chr_len)+'\t'+str(x) )
			if any(array<0): error ("array has reads less than 0")
			array = numpy.floor(array * scale_dict[filename])
			array = array + numpy.roll(array,-1)
	#		info("The total number of reads for %s after normalization is %s",filename,numpy.sum(array))
			if idx ==0:
				chr_array = array
			else:
				chr_array = numpy.column_stack((chr_array,array))
		data_by_window_dict[chr] = chr_array
	readData.reads_dict = data_by_window_dict

