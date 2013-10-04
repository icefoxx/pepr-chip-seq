from optparse import OptionParser
from classDef import ReadData


def optParser(argv):
	parser = OptionParser()
	parser.add_option("-c","--chip",action="store",
			 type="string",dest="chip",default="",
			 help="chip file names separated by comma",metavar="CHIP")
	parser.add_option("-i","--input",action="store",
			 type="string",dest="control",default="",
			 help="input/control file names separated by comma",metavar="INPUT")
	parser.add_option("-f","--fileformat",action="store",
			 type="string",dest="fileFormat",
			 help="file formats: bed,eland_multi,eland_extended,sam,bowtie...", metavar="FORMAT")
	parser.add_option("-g","--genome",action="store",
			 type="string",dest="species",
			 help="mm9,hg19,hg18,...",metavar="SPECIES")
	parser.add_option("-s","--shiftsize",action="store",
			 type="string", dest="shiftSize",default="-1",
			 help="Half the fragment size. If not provided, we'll estimate it",
			 metavar="SHIFTSIZE")
	parser.add_option("-w","--windowsize",action="store",
			 type="int",dest="windowSize",default=-1,
			 help="If not provided, we'll will estimate it",
			 metavar="WINDOWSIZE")
	parser.add_option("-n","--name",action = "store",
			 type="string",dest="name",default = "NA",
			 help = "the experimental name",
			 metavar="NAME")
	parser.add_option("-r","--remove_duplicate",action ="store_true",
			 dest = "remove_redundant",default=False,
			 help="Whether to remove duplicated reads")
	parser.add_option("--threshold",action ="store",
			 type='float',dest="threshold",default=1e-5)
	parser.add_option("--peaktype",action="store",
			 type="string",dest="peaktype",default="broad",
			 help="sharp or broad.")
	parser.add_option("--filter_short_fragment",action="store_true",
			 dest="remove_short_fragment",default=False,
			 help = "Whether to remove peaks that have estimated fragment size less than readLength+5. Only available for SHARP peak type")
	parser.add_option("--narrow_peak_width",action="store_true",
			 dest ="narrow_peak_width",default=False,
			 help = "Whether to narrow peak width to contain the most enriched regions. Only available for SHARP peak type")
	(opt,args)=parser.parse_args(argv)

	return opt

def validateOpt(opt):
	if opt.chip =="":
		print " no test input detected. Try \"PePr -h\" for help"
		exit(1)
	if opt.control =="":
		print " no control input detected"
		exit(1)
	if opt.fileFormat not in ['bed','eland_multi', 'eland_extended', 'bowtie', 'sam']:
		raise Exception("file format input error")

	opt.peaktype = opt.peaktype.lower()
	if opt.peaktype not in ['sharp','broad']:
		raise Exception("please specify a peak type: sharp or broad. Typically, sharp works for TF better and broad for histone modifications.")

	## parse the chip and control filename 
	chip_filename_list=opt.chip.strip().split(',')
	control_filename_list=opt.control.strip().split(',')
	## initialize the data structure
	opt.read_data = ReadData(chip_filename_list,control_filename_list,opt.species)	
