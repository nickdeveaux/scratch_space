#!/usr/bin/python

# Author: Emily Mi

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
'GSM2056292.bed'
,'GSM2056293.bed'
,'GSM2056294.bed'
,'GSM2056295.bed'
,'GSM2056296.bed'
,'GSM2056297.bed'
,'GSM2056298.bed'
,'GSM2056299.bed'
,'GSM2056300.bed'
,'GSM2056301.bed'
,'GSM2056302.bed'
,'GSM2056303.bed'
,'GSM2056304.bed'
,'GSM2056305.bed'
,'GSM2056306.bed'
,'GSM2056307.bed'
,'GSM2056308.bed'
,'GSM2056309.bed'
,'GSM2056310.bed'
,'GSM2056311.bed'
,'GSM2056312.bed'
,'GSM2056313.bed'
,'GSM2056314.bed'
,'GSM2056315.bed'
,'GSM2056316.bed'
,'GSM2056317.bed'
,'GSM2056318.bed'
,'GSM2056319.bed'
,'GSM2056320.bed'
,'GSM2056321.bed'
,'GSM2056322.bed'
,'GSM2056323.bed'
,'GSM2056324.bed'
,'GSM2056325.bed'
,'GSM2056326.bed'
,'GSM2056327.bed'
,'GSM2056328.bed'
,'GSM2056329.bed'
,'GSM2056330.bed'
,'GSM2056331.bed'
,'GSM2056332.bed'
,'GSM2056333.bed'
,'GSM2056334.bed'
,'GSM2056335.bed'
,'GSM2056336.bed'
,'GSM2056337.bed')

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