
# initialize libraries.
import numpy
import logging
from classDef import Peak
from scipy.stats.distributions import poisson
from scipy.special import psi
from scipy.optimize import fsolve
from scipy.stats.distributions import norm
from operator import attrgetter

root_logger = logging.getLogger("")
debug = root_logger.debug
info = root_logger.info


def get_candidate_window(mean,actual,threshold):
	# get candidate windows which have a significant pvalue using a poisson model 
 	pvalue = 1-poisson.cdf(actual,mean)
	idx_list = numpy.where(pvalue[10:-10]<threshold)[0]+10
	return idx_list

def weighted_log_likelihood(v_hat, m, n, reads):
	equation = 0
	for idx in range(21):
		x = reads[idx,0:m]
		y = reads[idx,m:(m+n)]
		u_hat = numpy.mean(x)
		r_hat = numpy.mean(y)/u_hat

		log_likelihood = -(m+n)*psi(v_hat)+numpy.sum(psi(v_hat+x))+numpy.sum(psi(v_hat+y))+m*numpy.log(v_hat/(v_hat+numpy.mean(x)))+n*numpy.log(v_hat/(v_hat+numpy.mean(y)))

		equation = equation+log_likelihood  * (1- (abs(11-idx)/11)) #weight
	return equation



def cal_area_dispersion_factor(read,m,n, idx):
	# calculate an area dispersion factor for each window
	return fsolve(weighted_log_likelihood, 1, args=(m,n,read[(idx-10):(idx+11)]))[0]  #find the optimum dispersion parameter that maximize the likelihood. 

def cal_local_dispersion_factor(read,m,n, idx):
	return fsolve(weighted_log_likelihood2, 1, args=(m,n,read[idx]))

def weighted_log_likelihood2(v_hat, m, n, reads):
	x = reads[0:m]
	y = reads[m:(m+n)]
	log_likelihood = -(m+n)*psi(v_hat)+numpy.sum(psi(v_hat+x))+numpy.sum(psi(v_hat+y))+m*numpy.log(v_hat/(v_hat+numpy.mean(x)))+n*numpy.log(v_hat/(v_hat+numpy.mean(y)))

        return log_likelihood

def cal_FDR(peak_list,num_tests):
	
	peak_list = sorted(peak_list,key=attrgetter('pvalue'))
	#calculate BH q-values
	q_list = [item.pvalue*num_tests/(idx+1) for idx,item in enumerate(peak_list)]
	
	for i in range(len(q_list)):
		for j in range(i):
			if q_list[j] > q_list[i]:
				q_list[j]=q_list[i]
	for i in range(len(q_list)):
		peak_list[i] = Peak(peak_list[i].chr,peak_list[i].index,peak_list[i].pvalue,q_list[i])
	peak_list = sorted(peak_list, key=attrgetter('chr','index'))
	return peak_list

def negative_binomial(readData,peakfilename,peaktype,swap=False,remove_shift=False,narrow_peak=False,threshold=1e-5,windowsize=200):
	#the main function that test for significant windows. 

	read_dict = readData.reads_dict 
	strands_dict = readData.data_dict_by_strands
	chr_list = readData.chr_list
	chip_list = readData.chip_filename_list
	control_list = readData.control_filename_list
	num_tests = readData.genomeSize/windowsize
	peakfile = open(peakfilename,'w')

	#compute number of replicates
	chip_rep = len(chip_list)
	control_rep = len(control_list)
	start1 = 0
	end1 = start2 = chip_rep
	end2 = chip_rep+control_rep
	# if swap 
	if swap: 
		chip_rep,control_rep = control_rep,chip_rep
		chip_list,control_list = control_list,chip_list

	# initialize basic array structures
	#sig_index_dict = {}
	sig_peaks_list = []	

	for chr in chr_list:

		read_array = read_dict[chr]
		read_array[numpy.where(read_array ==0)] = 1 
		y_bar_array = numpy.mean(read_array[:,start1:end1],1)
		x_bar_array = numpy.mean(read_array[:,start2:end2],1)
		if swap: #swap the chip and control reads. 
			x_bar_array,y_bar_array = y_bar_array,x_bar_array
		# setting the minimum # of reads in each window to 1 so that they won't have arithmetic errors. 
	
		cand_index = get_candidate_window(x_bar_array,y_bar_array,threshold)  # define a window as candidate if its Poisson p-value is smaller than the threshold. 
		debug("there are %d candidate windows for %s",len(cand_index),chr)
		debug("begin estimating dispersion parameters")
		if not swap: 
			disp_list = numpy.array([cal_area_dispersion_factor(read_array, chip_rep, control_rep, idx) for idx in cand_index]) #for all candidate windows, calculate the dispersion parameters. 
		else: 
			disp_list = numpy.array([cal_area_dispersion_factor(read_array, control_rep, chip_rep, idx) for idx in cand_index])
		debug("finished estimating dispersion")
		cand_x_bar_array = x_bar_array[cand_index]
		cand_y_bar_array = y_bar_array[cand_index]
		gamma_array = cand_y_bar_array / cand_x_bar_array
		tau_hat_array = numpy.sqrt(cand_y_bar_array*((control_rep*cand_x_bar_array*(disp_list+cand_y_bar_array)) + (chip_rep*cand_y_bar_array*(disp_list+cand_x_bar_array))) / (chip_rep*control_rep*disp_list*(cand_x_bar_array**3)))

		gamma_hat = 1.0 #Null hypothesis
		z_score_array = (numpy.log(gamma_array)-numpy.log(gamma_hat))*gamma_array/tau_hat_array

		pval_array = norm.cdf(-z_score_array)

		pois_pvalue = 1-poisson.cdf(cand_y_bar_array,cand_x_bar_array)

		# write the poisson output files. for debugging purposes. 
#                file2 = open(chr+".poisson",'w')
#                for i,idx in enumerate(sig_index):
#                        file2.write(str(idx *windowsize/2)+'\t'+str(idx*windowsize/2+windowsize)+'\t')
#                        read = read_array[idx]
#                        for item in read:
#                                file2.write(str(item)+'\t')
#                        file2.write(str(pois_pvalue[i])+'\t'+str(disp_list[i])+'\t'+str(pval_array[i])+'\n')


		test_index = numpy.where(pval_array<threshold)  # record the indices of the windows that have p-value smaller than the threshold. 
		test_index = test_index[0]
		sig_index = cand_index[test_index]
		sig_pval = pval_array[test_index]
		sig_disp = disp_list[test_index]
		for i,a in enumerate(test_index):
			sig_peaks_list.append(Peak(chr,sig_index[i],sig_pval[i],0))

	#calculate the BH FDR. 
	debug("begin estimating fdr")
	sig_peaks_list = cal_FDR(sig_peaks_list,num_tests)
	debug("finished estimating fdr")

	for chr in chr_list:
		sig_peak_list_by_chr = [item for item in sig_peaks_list if item.chr==chr]
		sig_index = [item.index for item in sig_peak_list_by_chr]
		sig_pval = [item.pvalue for item in sig_peak_list_by_chr]
		sig_qval = [item.qvalue for item in sig_peak_list_by_chr]	 
		sig_start,sig_end,sig_pval,sig_qval = merge_sig_window(sig_index,sig_pval,sig_qval,peaktype)
			
		if peaktype =="sharp":
			for idx in xrange(len(sig_start)):
				decision = shift_size_per_peak(strands_dict, chip_list, chr, sig_start[idx]*windowsize/2, sig_end[idx]*windowsize/2+windowsize,readData.shiftSize,readData.readLength,narrow_peak)
				if (decision==-1):
					continue   #skip this iteration; means that this peak is skipped and not reported in the final peak list. 
				start,end,shift_list,foldc = decision
				if ((validate_shift_size(shift_list,readData.readLength)==0) & remove_shift):
					continue  #skip this iteration; means that this peak is skipped and not reported in the final peak list.
				peakfile.write(chr+"\t"+str(start)+'\t'+str(end)+'\t'+str(end-start)+'\t')

				for shiftsize in shift_list:
					peakfile.write(str(shiftsize)+'\t')
				#peakfile.write(str(foldc)+'\t') #report the ratio of the forward strand reads over the reverse strand reads. #will not report in this version.
				peakfile.write(str(sig_pval[idx])+'\t'+str(sig_qval[idx])+'\n')

		elif peaktype == "broad":
			for idx in xrange(len(sig_start)):
				peakfile.write(chr+"\t"+str(sig_start[idx]*windowsize/2)+'\t'+str(sig_end[idx]*windowsize/2+windowsize)+'\t'+str((sig_end[idx]-sig_start[idx])*windowsize/2+windowsize)+'\t')
				peakfile.write(str(sig_pval[idx])+'\t'+str(sig_qval[idx])+'\n')
	
	return

def merge_sig_window(index_list,pval_list,qval_list,peaktype):
	# merge significant windows that are nearby. 
	if peaktype=="sharp":
		MIN_WINDOW = 1 #mininal number of windows required. 
		MAX_WINDOW = 2 #the maximal gap (max-1) of significant windows allowed to merge. 
	elif peaktype =="broad":
		MIN_WINDOW = 2 
		MAX_WINDOW = 5
	sig_peak_start = []
	sig_peak_end = []
	sig_pval = []
	sig_qval = []
	for idx,pos in enumerate(index_list):
		if idx == 0:
			sig_peak_start = []
			sig_peak_end = []
			start = pos
			pre_start = pos
			peak_pval_list = [pval_list[idx]]
			peak_qval_list = [qval_list[idx]]
		else: 
			if pos - pre_start <= MAX_WINDOW:
				pre_start = pos
				peak_pval_list.append(pval_list[idx])
				peak_qval_list.append(qval_list[idx])

			elif pos-pre_start > MAX_WINDOW:
				end = pre_start
				if (end-start>=MIN_WINDOW-1):
					sig_peak_start.append(start)
					sig_peak_end.append(pre_start)
					sig_pval.append(min(peak_pval_list))
					sig_qval.append(min(peak_qval_list))
				start = pos
				pre_start = pos
				peak_pval_list = [pval_list[idx]]
				peak_qval_list = [qval_list[idx]]

	if (pre_start-start>MIN_WINDOW-1):
		sig_peak_start.append(start)
		sig_peak_end.append(pre_start)
		sig_pval.append(min(peak_pval_list))
		sig_qval.append(min(peak_qval_list))

	return sig_peak_start,sig_peak_end,sig_pval,sig_qval

def shift_size_per_peak(strands_dict, chip_list, chr, start, end, shiftSize,readLength,narrow_peak):
	shift_list = []
	start_list = []
	end_list = []

	for chip in chip_list:
		forward = strands_dict[chr][chip]['f']
		reverse = strands_dict[chr][chip]['r']

		forward_read = forward[numpy.where( (forward > start-shiftSize[chip]) & (forward < end-shiftSize[chip]) )]
		reverse_read = reverse[numpy.where( (reverse > start+shiftSize[chip]) & (reverse < end+shiftSize[chip]) )]
		if len(reverse_read) ==0: 
			fold = 100  #if there is no reads in the reverse strand, define the fold (forward over reverse) to be 100 to avoid arithmetic errors. 
		else:
			fold = float(len(forward_read))/len(reverse_read)
		# estimate the shift size
		forward = numpy.zeros(end-start)
		reverse = numpy.zeros(end-start)
		for read in forward_read:
			forward[(read-start):(read-start+readLength)] +=1
		for read in reverse_read:
			reverse[(read-start-readLength):(read-start)] +=1
		shade_max = 0
		shift_max = 0
		for shift in xrange(-readLength,250-readLength): #iterate from 0 to 250 bp, find the optimum shift that produce the most overlap between the two strands. 
			forward_temp = numpy.roll(forward,shift)
			shade = numpy.sum(numpy.min(numpy.vstack([forward_temp,reverse]),0 ) )
			if shade > shade_max:
				shift_max = shift
				shade_max = shade
		
		shift_list.append(shift_max+readLength)
		if narrow_peak:
			forward_read.sort()
			reverse_read.sort()
			forward_len = len(forward_read)
			reverse_len = len(reverse_read)
			if forward_len>0:
				start_list.append(forward_read[numpy.floor(forward_len*0.2)]) #refine the left boundary of the peak.  
			if reverse_len>0: 
				end_list.append(reverse_read[-numpy.floor(reverse_len*0.2)])  #refine the right boundary of the peak. 
	if narrow_peak:
		if (len(start_list) ==0) | (len(end_list)==0):
			return -1
		else:
			return min(start_list),max(end_list),shift_list,fold
	else:
		return start,end,shift_list,fold


def validate_shift_size(list,readLength):
	#if the shift size is smaller than the readlength+5, then the peak will be removed. 
	list_len = len(list)
	short_len =0
	for shift in list:
		if shift <=readLength+5:
			short_len = short_len+1

	if short_len > list_len/2:
		return 0
	else:
		return 1

