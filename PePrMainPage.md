The latest PePr manual is updated [here](http://ones.ccmb.med.umich.edu/wiki/PePr/).



# Introduction #

ChIP-seq is a standard method to identify genome-wide protein-DNA interaction or histone modification sites. A substantial number of peak-calling programs exist to automatically identify binding sites from ChIP-Seq data; PePr differs in that it takes into account variance among biological replicates. Even with transcription factor experiments, accounting for variation has the potential to identify binding sites more consistently functionally relevant in the context under study.

Here we offer a novel ChIP-Seq Peak-calling and Prioritization pipeline (PePr) that uses a sliding window approach and models read counts across replicates and between groups with a negative binomial distribution. PePr empirically estimates the optimal shift/fragment size and sliding window width, and estimates dispersion from the local genomic area. Regions with less variability across replicates are ranked more favorably than regions with greater variability. Optional post-processing steps are also made available to filter out peaks not exhibiting the expected shift size and/or to narrow the width of peaks.


---


# Installation #
  1. Make sure your python version is 2.6 or another that is compatible with the numpy and scipy packages.
  1. **(Highly recommended)** If you have installed pip in your system, simply type  `"pip install PePr"`  or  `"pip install PePr --user"`  (if you don't have administrator privilege).
  1. If you don't want to use pip, you have to manually download and install dependencies (numpy, scipy, and pysam packages). And then download **[PePr](https://pypi.python.org/pypi/PePr)**. Type `"tar zxvf PePr-[version].tar.gz"` to uncompress the file. Then go into the PePr directory and type `"python setup.py install --prefix=~/.local/"` to install it in your home directory
  1. Go into the python package directory (Depending on where you install it, it may be $HOME/.local/lib/pythonX.Y/site-packages/PePr/) and type “python ./PePr.py –h” to verify that the help messages pop up. In addition, go into the PePr/data directory to run the following test command to see if they produce the expected results.
```
python ../PePr.py -i input_rep1.bed,input_rep2.bed -c chip_rep1.bed,chip_rep2.bed -f bed -s 48,46 -w 180 -n my_test_run
```
  1. If they do, you’re ready to go!
  1. To make life easier, you can create a short-cut to PePr using alias command in linux: `alias pepr="python thePathToPePr/PePr/PePr.py"`. And put that command in your $HOME/.bashrc file.

---

# How to run PePr #
There are some minor differences between running a peak-calling pipeline and a differential binding analysis. For peak-calling, only three parameter fields are currently required, `-c`, `-i` and `-f`, which are ChIP files, input files and file format respectively. For differential binding analysis, three additional parameter fields are required, which are `--chip2`, `--input2` and `--diff`. `--chip2` and `--input2` specify the ChIP and input files of the other condition, while `--diff` tells PePr to run in "differential binding" mode. These are minimum inputs for PePr to begin analysis.
Additionally, you are encouraged to specify the experiment name (otherwise a default "NA" will be used). Attn: if there are multiple runs of PePr with the same name in the same folder, older results will be replaced. It is also important to tell PePr what kind of shape you expect for your data, either sharp (like Transcription Factors, H3K4me3) or broad (like H3K27me3). For all datasets, it is strongly recommended to leave `--remove_artefacts` on. It will remove PCR-duplicated peaks from the peak list and report them in a different file (so you will not lose that peak if PePr made a mistake(oops! sorry!)). If you're not familiar with this kind of false positive we're talking about, please check out the example [here](ExampleNRSF#4._Removing_artefacts.md)

  1. For transcription factors, we recommend the following options:
```
python path/PePr.py -i inputFiles -c chipFiles -n ExperimentName -f format --peaktype=sharp --remove_artefacts
```
  1. For histone modifications that are expected to cover large genomic regions, we recommend the following options:
```
python path/PePr.py -i inputFiles -c chipFiles -n ExperimentName -f format --peaktype=broad  --remove_artefacts
```
  1. For differential binding analysis of transcription factors:
```
python path/PePr.py -i inputFiles_group1 --input2 inputFiles_group2 -c chipFiles_group1 --chip2 chipFiles_group2 -n ExperrrimentName -f format --peaktype=sharp --diff --remove_artefacts
```
  1. For differential binding analysis of transcription factors without input samples:
```
python path/PePr.py -c chipFiles_group1 --chip2 chipFiles_group2 -n ExperrrimentName -f format --peaktype=sharp --diff --remove_artefacts
```
  1. Example of exact command for a transcription factor with two input files and two ChIP files:
```
python path/PePr.py –i input-1.bed,input-2.bed –c chip-1.bed,chip-2.bed –n myTF –f bed --peaktype=sharp –-remove_artefacts 
```

---

# Input format #
**bed**, **sam**, **bam**, eland\_multi, eland\_exetnded and bowtie format are currently supported.

Note: In order for **bam** files to load, user will need to sort and index them first. You may use the following commands:
```
samtools sort your_file.bam your_file.sorted 
samtools index your_file.sorted.bam
```

---

# Parameters #
**-h, --help**
> Show this help message and exit.
**-i INPUT1**
> Group 1 input files. **REQUIRED.** Multiple file names are
> separated by comma,e.g. input1.bed,input2.bed
**-c CHIP1**
> ChIP-Seq chip group 1 files. **REQUIRED.** Multiple file names are
> separated by comma, e.g. chip1.bed,chip2.bed
**--input2 INPUT2**
> Group 2 input files. Used in differential binding analysis.
**--chip2 CHIP2**
> Group 2 ChIP-Seq files. Used in differential binding analysis.
**-n NAME**
> Experiment name. DEFAULT: "NA"
**-f FORMAT**
> File format. **REQUIRED.** (currently support bed,
> eland\_multi, eland\_exetnded, bowtie, sam, bam)
**-s SHIFTSIZE**
> Half fragment size. The number of bases to shift
> forward and reverse strand reads toward each other.
> If not specified by user, PePr will empirically
> estimate this number from the data. A separate shift
> size will be estimated for each ChIP file.
**-w WINDOWSIZE**
> Window size. If not specified by user, PePr will
> estimate this from the data.
**--diff**
> Whether to perform differential binding analysis.
> If so, please provide chip2 and input2 samples as well.
**--threshold THRESHOLD**
> Significance cutoff for p-value DEFAULT:1e-5
**-r, –-remove\_duplicate**
> If this option is specified, redundant reads greater
> than the expected maximum (determined by a binomial
> test) at each single location will be removed.
**--peaktype PEAKTYPE**
> Is peak shape sharp or broad?(sharp, broad)
**--remove\_artefacts**
> Remove peaks that are caused by PCR duplicates.
**--narrow\_peak\_width**
> Narrow peak width to contain the most enriched region.
> Only available for SHARP peak type.

---

# Output Files #

  * **NAMEPePr\_peaks.bed:** A tab-delimited file containing chromosomal position of the peak, peak width, p-value and Benjamini-Hochberg FDR.
  * When `diff` option is specified, you will see NAMEPePr\_up\_peaks.bed and NAMEPePr\_down\_peaks.bed, which shows peaks in chip1 and chip2 samples respectively.
> > Example columns for peaks:

| **chrom** | **start** | **end** | **peak width** | **p-value** | **q-value** |
|:----------|:----------|:--------|:---------------|:------------|:------------|
|chr1       |837200     |837500   |300             |1.29969869177e-18|3.12684869127e-15|
|chr1       |960800     |961100   |300             |1.58064616217e-09|1.93155073885e-06|
|chr1       |1224000    |1225400  |1400            |3.45692322029e-86|2.9886376305e-80|

  * **NAME\_PePr\_filtered\_peaks.txt:** A tab-delimited file showing the peaks filtered by PePr. Only present when --remove\_artefacts option is specified.

  * **NAME\_PePr\_parameters.txt:** A file containing the parameters to reproduce the results.
  * **Date-debug.log:** This file contains the detailed information about the running status. Useful debugging information contains: Total reads per chrom, maximum number of reads at a single position by binomial test, etc. (miscellaneous advanced details).

---

# Contact #
For more information about our work please contact:
> > Maureen Sartor, Ph.D. - sartorma@med.umich.edu [Lab website](http://sartorlab.ccmb.med.umich.edu/)
General question and report bugs please contact:
> > Yanxiao Zhang - yanxiazh@umich.edu [Personal website](http://www-personal.umich.edu/~yanxiazh/)

If you like this, give us a +1:) 