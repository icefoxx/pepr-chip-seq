# Introduction #

This page is intended as a comprehensive example to illustrate the essential processing steps of PePr (See [Methods section](PePrMethods.md) for details) with the **peak-calling** setting.

# Data description #
Neuron-Restrictive Silencer Factor (NRSF) is a ENCODE TF dataset. It has two replicates each for both ChIP and input sample groups. The following links lead to the data download:
  * [chip rep1](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeHudsonalphaChipSeq/wgEncodeHudsonalphaChipSeqAlignmentsRep1K562Nrsf.tagAlign.gz)
  * [chip rep2](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeHudsonalphaChipSeq/wgEncodeHudsonalphaChipSeqAlignmentsRep2K562Nrsf.tagAlign.gz)
  * [input rep1](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeHudsonalphaChipSeq/wgEncodeHudsonalphaChipSeqAlignmentsRep1K562Control.tagAlign.gz)
  * [input2 rep2](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeHudsonalphaChipSeq/wgEncodeHudsonalphaChipSeqAlignmentsRep2K562Control.tagAlign.gz)

The files are in tagAlign format, which is very similar to BED format and can be processed in the same way as BED format by PePr. The total # of aligned reads are summarized below:
| **Sample** | **Number of Mapped Reads (million)** |
|:-----------|:-------------------------------------|
| chip-rep1  | 16.1                                 |
| chip-rep2  | 26.6                                 |
| input-rep1 | 16.3                                 |
| input-rep2 | 14.3                                 |

Note that ChIP replicate 1 has many more reads as replicate 2.

# PePr analysis on this data #
In the following subsections, we will show the command we executed on this data, the pre- and post-processing steps.

## commands: ##
Assuming we have renamed the ChIP samples and input samples accordingly.
```
python pathToPePr/PePr.py -i input-1.tagAlign,input-2.tagAlign -c chip-1.tagAlign,chip-2.tagAlign
-f bed --peaktype=sharp --remove_artefacts
```
## processing steps ##
### 1. Shift/fragment size estimates ###
PePr estimated the shift sizes to be 46/48 for ChIP rep1/rep2. The following figures show the distribution of shift/fragment sizes estimated from the peak regions of these two replicates.

![http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_shift_size_distribution-rep1.png](http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_shift_size_distribution-rep1.png)
![http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_shift_size_distribution-rep2.png](http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_shift_size_distribution-rep2.png)

The small bump on the left that is around the read length is the "phantom peak", which represents the peaks caused by high level of PCR-duplicates.

### 2. Window size estimates ###
PePr's optimal window size estimate for this data is 200bp. Below shows a representative peak region (Using IGV for visualization) that have a width similar to PePr's estimtes.

![http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_window_size_estimate.png](http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_window_size_estimate.png)

### 3. Normalization ###
As there are unequal number of reads in every sample, we need to normalize for the differences in sequencing depth and ChIP efficiencies. PePr used the mean of all ChIP replicates as the reference sample, and normalize every other sample against it (For details please see [Methods section](Methods.md).
#### (1) Normalization of input ####
For the normalization of input sample, PePr uses the NCIS method. First the whole genome is split into 1KB windows. For each read threshold r (from 1 to 400), we find the windows having reads less than r, and plot the stats below (%genome these windows cover, %reads these windows have compared to whole libaray, and ratio of reference over input):
![http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_input_normalization.png](http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_input_normalization.png)
In this figure, the red points show the PePr's estimate of the normalization constants (on the third subplots), which pretty much stabilizes after r>200. The blue dotted lines show the ratio of total library reads.
#### (2) Normalization of ChIP ####
The ChIP samples are normalized against the above-mentioned reference sample as well. First the whole genome are split into 1KB windows. Then the windows are sorted by the total number of reads in the window. We don't know which windows are enriched and which are simply background (some windows are a mixture of both); however, we can assume that windows with more reads are more likely to be enriched than background. The goal in this step is to estimated the signal ratio (or ChIP efficiency ratio) of the enriched regions. So we estimate the normalization constants for the top N windows where N ranges from 1000 to 50000, to see how the ratio changes across different number of windows we include.
![http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_chip_normalization.png](http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/NRSF_chip_normalization.png)
The figure above shows the ratio versus N in this data. The blue dotted line is the library ratio. We can see that the normalization constant first depart from the library ratio and approaches a maximum/minimum when N increases (where all windows are likely enriched region), and then converge back to library ratio as N further increases (beyond the point where the windows are likely background).
### 4. Removing artefacts ###
PePr removes the peaks that are highly likely PCR-dupliation artefacts.
The below is an example of the peaks that are removed by PePr.
![http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/chr4-49349200-49349500.png](http://www-personal.umich.edu/~yanxiazh/PePr/exampleNRSF/chr4-49349200-49349500.png)