import GenomeData
import logging
import misc

rootlogger = logging.getLogger("")
info = rootlogger.info
debug = rootlogger.debug

class ReadData:
	''' A data structure for the read data'''
	def __init__(self,chip,control,sp):
		self.data_dict_by_strands = {}
		self.data_dict = {}
		self.reads_dict = {}
		self.chip_filename_list = chip
		self.control_filename_list = control
		self.filename_list =chip+control
		self.chr_list, self.chr_length_dict = GenomeData.species_match(sp)
		self.read_total_per_file = {}
		self.genomeSize = 0
		self.shiftSize = {}
		self.readLength = 0
	def cal_genomeSize(self):  #calculate the genome size.
		for chr in self.chr_list:
			self.genomeSize += self.chr_length_dict[chr]
		debug("The genome size is %d",self.genomeSize)
	def cal_read_total(self): # calculate the total number of reads for each sample.	
		for  chr in self.chr_list:
			for filename in self.filename_list:
				for strand in ["f","r"]:
					try: self.read_total_per_file[filename] += len(self.data_dict_by_strands[chr][filename][strand])
					except KeyError: self.read_total_per_file[filename] = len(self.data_dict_by_strands[chr][filename][strand])

	def remove_redundant_reads(self):
		#remove additional redundant reads that are not warranted by a binomial test. 

		def rm_max(list,max):
			#remove reads that are more than the maximum specified. 
			list_new = []
			pre = 0
			count = 0

			for x in list:
				if x == pre: 
					count += 1
				else:
					pre = x
					count = 1
				if count > max:
					continue
				list_new.append(x)
			return list_new
					

		for filename in self.filename_list:
			max = misc.binomial(self.read_total_per_file[filename],1.0/self.genomeSize) #calculate maximum duplicates for each genomic location.
			info (" The maximum allowed reads at one position for %s is %d",filename, max)
			total_before = self.read_total_per_file[filename]
			total_after = 0
			for chr in self.chr_list:
				for strand in ["f","r"]:
					reads = self.data_dict_by_strands[chr][filename][strand]
					reads.sort()
					reads = rm_max(reads,max)
					total_after += len(reads)
					self.data_dict_by_strands[chr][filename][strand] = reads	
			self.read_total_per_file[filename] = total_after
			debug( " Total # of reads before removing redundancies is %d",total_before)
			debug( " Total # of reads after removing redundancies is %d",total_after)
			info( "The percentage of redundant reads for %s is %f ",filename,1-float(total_after)/total_before)
			
	def refine_chr_list(self): #delete the chromosomes that having 0 reads. This is just for testing purposes when only one or few chromosome was used. 
		chr_list = []
		for x in self.chr_list:
			if self.data_dict_by_strands[x][self.filename_list[0]]['f'] == []:
				continue
			chr_list.append(x)
		self.chr_list=chr_list


class Peak:
	"data structure that contains the significant peaks" 

	def __init__(self,chr,index,pvalue,qvalue):
		self.chr = chr
		self.index = index
		self.pvalue = pvalue
		self.qvalue = qvalue
	def __repr__(self):
		return repr((self.chr,self.index,self.pvalue,self.qvalue))

