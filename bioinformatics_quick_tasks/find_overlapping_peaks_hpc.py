#!/usr/bin/python

Usage = """
find_overlapping_peaks.py
For a set of .bed files containing peaks, this script will make a similarity matrix, where
each pairwise similarity is the % overlap in peaks 
Output is an asymmetric matrix of % overlaps -- asymmetric because if multiple peaks in sample
A map to a single peak in map B, these will be counted a single time in entry (B,A) and multiple 
times in entry (A,B) of the matrix
USAGE:
INPUTS:
	OUTPUTfolder -- contains function output (one matrix per TF, as a text file containing
		accessibility quantification centered on predicted motif sites in accesible regions)
"""	

import HTSeq
import numpy
import pysam
import sys
import itertools
import re
from scipy.stats import hypergeom
import time
import math

beddirectory = '/mnt/ceph/users/ndeveaux/ILC_NERD_2017/peakdeck/peakdeck_75bps_10kb/'

bedfiles = (
'GSM2056292_sorted_and_merged.bed'
,'GSM2056293_sorted_and_merged.bed'
,'GSM2056294_sorted_and_merged.bed'
,'GSM2056295_sorted_and_merged.bed'
,'GSM2056296_sorted_and_merged.bed'
,'GSM2056297_sorted_and_merged.bed'
,'GSM2056298_sorted_and_merged.bed'
,'GSM2056299_sorted_and_merged.bed'
,'GSM2056300_sorted_and_merged.bed'
,'GSM2056301_sorted_and_merged.bed'
,'GSM2056302_sorted_and_merged.bed'
,'GSM2056303_sorted_and_merged.bed'
,'GSM2056304_sorted_and_merged.bed'
,'GSM2056305_sorted_and_merged.bed'
,'GSM2056306_sorted_and_merged.bed'
,'GSM2056307_sorted_and_merged.bed'
,'GSM2056308_sorted_and_merged.bed'
,'GSM2056309_sorted_and_merged.bed'
,'GSM2056310_sorted_and_merged.bed'
,'GSM2056311_sorted_and_merged.bed'
,'GSM2056312_sorted_and_merged.bed'
,'GSM2056313_sorted_and_merged.bed'
,'GSM2056314_sorted_and_merged.bed'
,'GSM2056315_sorted_and_merged.bed'
,'GSM2056316_sorted_and_merged.bed'
,'GSM2056317_sorted_and_merged.bed'
,'GSM2056318_sorted_and_merged.bed'
,'GSM2056319_sorted_and_merged.bed'
,'GSM2056320_sorted_and_merged.bed'
,'GSM2056321_sorted_and_merged.bed'
,'GSM2056322_sorted_and_merged.bed'
,'GSM2056323_sorted_and_merged.bed'
,'GSM2056324_sorted_and_merged.bed'
,'GSM2056325_sorted_and_merged.bed'
,'GSM2056326_sorted_and_merged.bed'
,'GSM2056327_sorted_and_merged.bed'
,'GSM2056328_sorted_and_merged.bed'
,'GSM2056329_sorted_and_merged.bed'
,'GSM2056330_sorted_and_merged.bed'
,'GSM2056331_sorted_and_merged.bed'
,'GSM2056332_sorted_and_merged.bed'
,'GSM2056333_sorted_and_merged.bed'
,'GSM2056334_sorted_and_merged.bed'
,'GSM2056335_sorted_and_merged.bed'
,'GSM2056336_sorted_and_merged.bed'
,'GSM2056337_sorted_and_merged.bed')

maxpeaks = 300000
overlapmin = .75

OUTPUTfolder = '/mnt/ceph/users/ndeveaux/ILC_NERD_2017/Atac_Stats/'

## END INPUTS

output = open(OUTPUTfolder + 'peak_overlaps_top' + str(maxpeaks) + '_ov' + \
	str(overlapmin).strip('.') + '.txt','w')
covoutput = open(OUTPUTfolder + 'basecoverage_top' + str(maxpeaks) + '_ov' + \
	str(overlapmin).strip('.') + '.txt','w')
covoutput.write('Sample\tBasesInPeaks\n')

# get genomic arrays
peakregions = dict()
peakcutoffs = dict()
samplenames = dict()
maxpeaklength = dict()
for count in range(0,len(bedfiles)):
	samp = bedfiles[count]
	# get sample name
	samplenames[samp] = samp.strip('noM_peaks.bed')
	output.write('\t' + samp.strip('noM_peaks.bed'))
	# get peak locations		
	peakregions[samp] = HTSeq.GenomicArray("auto",stranded=False,typecode='d')
	maxpeaklength[samp] = 0
	scoreslist = list()
	peakfile = HTSeq.BED_Reader(beddirectory + samp)
	totbases = 0
	for peak in peakfile:
		peakregions[samp][peak.iv] = peak.score
		scoreslist.append(peak.score)
		peaklength = peak.iv.end - peak.iv.start + 1
		totbases += peaklength
		maxpeaklength[samp] = max(maxpeaklength[samp],peaklength)
	# find score cutoff for this particular library
	sortedscores = sorted(scoreslist)
	scorecutind = len(sortedscores) - maxpeaks - 1 # -1 for python indexing 0
	peakcutoffs[samp] = sortedscores[scorecutind]
	covoutput.write(samp.strip('noM_peaks.bed') + '\t' + str(totbases) + '\n')
covoutput.close()	
output.write('\n')

# now calculate overlaps
for count1 in range(0,len(bedfiles)):
	samp1 = bedfiles[count1]
	samp1cutoff = peakcutoffs[samp1]
	output.write(samplenames[samp1]) # + '\t0'*(count1+1))
	for count2 in range(0,len(bedfiles)):
		samp2 = bedfiles[count2]
		samp2cutoff = peakcutoffs[samp2]
		halfmaxshoulder = int(float(maxpeaklength[samp2])/2)
		curroverlaps = 0
		for peak1,score1 in peakregions[samp1].steps():
			if score1 >= samp1cutoff:	# is in top X significant peaks
				foundoverlap = 0
				# see if the peak also appears in the other sample
				for peak2,score2 in peakregions[samp2][HTSeq.GenomicInterval(peak1.chrom,\
					max(0,peak1.start-halfmaxshoulder),peak1.end+halfmaxshoulder,'.')].steps():
					if (score2 >= samp2cutoff and foundoverlap ==0) :
						if peak1.overlaps(peak2): # find out if they overlap 
							# find out how much they overlap
							inv1range = float(len(range(peak1.start,peak1.end+1)))
							inv2range = float(len(range(peak2.start,peak2.end+1)))
							overlap = float(len(range(max(peak1.start,peak2.start),\
								min(peak1.end,peak2.end))))
# 							print str(inv1range) + ' '+  str(inv2range) +' ' + str(overlap)
							if max(overlap/inv1range,overlap/inv2range) > overlapmin:
								curroverlaps += 1
								foundoverlap = 1 # we already found one overlapping peak
# 								print curroverlaps
		output.write('\t' + str(float(curroverlaps)/maxpeaks))
	output.write('\n')	
		
		
output.close()																