from optparse import OptionParser
from classDef import ReadData
from classDef import Parameters


def opt_parser(argv):
    parser = OptionParser()
    parser.add_option(
            "-p", "--parameter-file", action="store", type="string",
            dest="parameter", default="", 
            help="provide a file that contain the parameters",
            metavar="PARAMETER-FILE")
    parser.add_option(
            "-c", "--chip1", action="store", type="string", dest="chip1",
            default="", help="chip1 file names separated by comma",
            metavar="CHIP1")
    parser.add_option(
            "-i", "--input1", action="store",
            type="string", dest="input1", default="",
            help="input1 file names separated by comma",
            metavar="INPUT1")
    parser.add_option(
            "--chip2", action="store", type="string", dest="chip2",
            default="", help="chip2 file names separated by comma",
            metavar="CHIP2")
    parser.add_option(
            "--input2", action="store", type="string", dest="input2",
            default="", help="input2 file names separated by comma",
            metavar="INPUT2")
    parser.add_option(
            "-f", "--file-format", action="store", type="string", 
            dest="file_format",
            help="bed, eland_multi, eland_extended, sam, bowtie...",
            metavar="FORMAT")
    parser.add_option(
            "-s", "--shiftsize", action="store",
            type="string", dest="shift_size", default="-1",
            help="Half the fragment size.", metavar="SHIFTSIZE")
    parser.add_option(
            "-w", "--windowsize", action="store",
            type="int", dest="window_size", default=-1,
            help="Window sizes",
            metavar="WINDOWSIZE")
    parser.add_option(
            "--diff", action="store_true",
            dest = "difftest", default=False,
            help="Whether to perform differential binding analysis")
    parser.add_option(
            "-n", "--name", action = "store",
            type="string", dest="name", default = "NA",
            help = "the experimental name",
            metavar="NAME")
    parser.add_option(
            "-r", "--remove_duplicate", action ="store_true",
            dest = "remove_redundant", default=False,
            help="Whether to remove duplicated reads")
    parser.add_option(
            "--threshold", action ="store",
            type='float', dest="threshold", default=1e-5)
    parser.add_option(
            "--peaktype", action="store",
            type="string", dest="peaktype", default="broad",
            help="sharp or broad.")
    parser.add_option(
            "--remove_artefacts", action="store_true",
            dest="remove_artefacts", default=False,
            help = '''Whether to remove peaks that have the same shape 
            as in the input samples.''')
    parser.add_option(
            "--narrow_peak_width", action="store_true",
            dest ="narrow_peak_width", default=False,
            help = '''Whether to narrow peak width to contain the most 
            enriched regions. Only available for SHARP peak type''')
    parser.add_option(
            "--no_log", action="store_true",
            dest = "unsave_log", default=False,
            help = "Whether to disable saving the log files")
    (opt, args)=parser.parse_args(argv)

    return opt

def process_opt(opt):
    ''' validate the parameters that the user specified'''
    parameter = Parameters(opt)

    ## parse the chip and input filename 
    chip1_filename_list = parameter.chip1.strip().split(',')
    input1_filename_list = parameter.input1.strip().split(',')
    chip2_filename_list = parameter.chip2.strip().split(',')
    input2_filename_list = parameter.input2.strip().split(',')
    ## initialize the data structure
    read_data = ReadData(
        chip1_filename_list, input1_filename_list,
        chip2_filename_list, input2_filename_list,
        parameter.difftest 
        )
    #add shift size validations
    return parameter, read_data
