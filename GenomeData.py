##This file contains the genome information for most widely used organisms. Currently we support Mouse(mm9),Drosophila(dm3) and Human(hg18,hg19). Additional genomes could be easily integrated. 

mm9_chr_length_dict = {'chr1': 197195432,'chr2': 181748087,
		'chr3': 159599783,	'chr4': 155630120,
		'chr5': 152537259,	'chr6': 149517037,
		'chr7': 152524553,	'chr8': 131738871,
		'chr9': 124076172,	'chr10': 129993255,
		'chr11': 121843856,	'chr12': 121257530,
		'chr13': 120284312,	'chr14': 125194864,
		'chr15': 103494974,	'chr16': 98319150,
		'chr17': 95272651,	'chr18': 90772031,
		'chr19': 61342430,	'chrX': 166650296,
		'chrY': 15902555	}
dm3_chr_length_dict = {'chrUextra':29004656, 'chr3R': 27905053,
		'chr3L':24543557,	'chr2L': 23011544,
		'chrX': 22422827,	'chr2R': 21146708,
		'chrU': 10049037,	'chr2RHet': 3288761,
		'chr3LHet': 2555491,	'chr3RHet': 2517507,
		'chr4': 1351857,	'chr2LHet': 368872,
		'chrYHet': 347038,	'chrXHet': 204112,
		'chrM': 19517		}
hg18_chr_length_dict = {'chr1':  247249719, 'chr2': 242951149,
		'chr3': 199501827,	'chr4':    191273063,
		'chr5': 180857866,	'chr6':    170899992,
		'chr7': 158821424,	'chrX':    154913754,
		'chr8': 146274826,	'chr9':    140273252,
		'chr10': 135374737,	'chr11':   134452384,
		'chr12': 132349534,	'chr13':   114142980,
		'chr14': 106368585,	'chr15':   100338915,
		'chr16': 88827254,	'chr17':   78774742,
		'chr18': 76117153,	'chr20':   62435964,
		'chrY':  57772954,	'chr19':   63811651,
		'chr22': 49691432,	'chr21':   46944323 }
hg19_chr_length_dict = {'chr1':  249250621, 'chr2': 243199373,
		'chr3': 198022430,	'chr4':    191154276,
		'chr5': 180915260,	'chr6':    171115067,
		'chr7': 159138663,	'chrX':    155270560,
		'chr8': 146364022,	'chr9':    141213431,
		'chr10': 135534747,	'chr11':   135006516,
		'chr12': 133851895,	'chr13':   115169878,
		'chr14': 107349540,	'chr15':   102531392,
		'chr16': 90354753,	'chr17':   81195210,
		'chr18': 78077248,	'chr20':   63025520,
		'chrY':  59373566,	'chr19':   59128983,
		'chr22': 51304566,	'chr21':   48129895 }


mm9_chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
dm3_chr_list = ['chrUextra','chr3R','chr3L','chr2L','chrX','chr2R','chrU','chr2RHet','chr3LHet','chr3RHet','chr4','chr2LHet','chrYHet','chrXHet','chrM']
hg18_chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
hg19_chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

species_dict={'mm9':(mm9_chr_list,mm9_chr_length_dict),'dm3':(dm3_chr_list,dm3_chr_length_dict), 'hg18':(hg18_chr_list,hg18_chr_length_dict),'hg19':(hg19_chr_list,hg19_chr_length_dict)}

def species_match(arg):
	if arg not in species_dict:
		raise Exception("species not recognized")
	else:
		return species_dict[arg]

def main():
	chr_list,chr_length = species_match('mm9')
	print chr_list,chr_length
if __name__ == '__main__':
	main()
