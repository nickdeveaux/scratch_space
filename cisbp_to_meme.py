#!/usr/bin/python

# Author: Nick De Veaux
# This script uses the MIT licensed gscripts from UCSD's Yeo Lab for parsing pwms

# This script creates a meme formatted flat file database for motifs from cisbp, transfac and encode databases
# Some non-obvious decisions that were made: 
#   - Motifs can map to multiple TFs. If they do, the data field after the motif name has the comma seperated TF names
#   - Encode motifs are only removed as duplicates if their PWM has a .95 correlation with a PWM for the same TF
#       Unlike the encode 2013 paper, where the cutoff was .75 and offsets were used to search for duplicates of different lengths
#   - Motifs whose TF name is not in the complete human gene list have been removed
#       This was a hard decision, because it removes prefix-style TF symbols like "RAR", 
#        but there was no clean automatable way to map them,
#        and the number of TFs we would gain motifs for represented approx 1-2% of all TFs,
#        so the gains were considered marginal at the cost of potentially introducing incorrect TF to motif links. 

# Reference Data Files:
# TF class was used as the canonical source of human TFs (http://www.edgar-wingender.de/huTF_classification.html)
# combined with the GO annotations of DNA binding proteins (http://geneontology.org/page/download-annotations)
# CisBP 2016 (v1.02) motifs were downloaded from http://cisbp.ccbr.utoronto.ca/bulk.php
# Encode motifs came from http://compbio.mit.edu/encode-motifs/motifs.txt

import argparse
import os
from string import Template
import numpy as np
import pandas as pd
from gscripts import pwm

class MyTemplate(Template):
    delimiter = '&'

parser = argparse.ArgumentParser(description="""takes a cisbp or rbpdb formatted file and converts it to a meme / fimo formatted pwm""")
 
parser.add_argument("--pwm_folder", "-p", help="cisbp/rbpdb pwm folder", required=True)
parser.add_argument("--out_file", "-o", help="output file", required=True)
parser.add_argument("--information_file", "-i", help="information file", required=True)
parser.add_argument("--debug", help="verbose debug mode", action='store_true')
parser.add_argument("--mouse", help="run for mouse", action='store_true')

args = parser.parse_args()

transfac_file = "/Users/ndeveaux/Data/motifs/transfac_output.meme"
if args.mouse:
    transfac_file = "/Users/ndeveaux/Data/motifs/mouse_transfac_output.meme"

encode_file = "/Users/ndeveaux/Data/motifs/Encode_Site/motifs_manually_touched_up.txt"
gene_list_file = '/Users/ndeveaux/Data/human_reference/genes.gtf'

with open(os.path.join(pwm.pwm_dir(), "cisbp_header.txt")) as input:
    header = input.read()

output = header

with open(os.path.join(pwm.pwm_dir(), "cisbp_template.txt")) as input:
    template = MyTemplate(input.read())

if args.information_file:
    tf_information = pd.read_csv(args.information_file, sep="\t", header=0)

pwm_to_tf_names = {}
motif_matrices = {}

# Due to licensing issues, the TRANSFAC motifs (that we have a license for!) are not able to be downloaded from CisBP
# The solution is to save their IDs to a file that will be based in via the -id flag to Meme Suite's transfac2meme tool
missing_motifs = {}

# Iterate through the directory of CisBP files, building the required data structures
for pwm_file in os.listdir(args.pwm_folder):
    try:
        pwm_name = ".".join(pwm_file.split(".")[:-1])
        if pwm_name in tf_information['Motif_ID'].tolist():
            if args.debug:
                print pwm_name
            pwm_to_tf_names[pwm_name] = []
            potential_transfac_id = None
            motif_rows = tf_information[tf_information['Motif_ID'] == pwm_name]
            for idx, row in motif_rows.iterrows():
                tf_name = row['TF_Name']
                pwm_to_tf_names[pwm_name].append(tf_name)
                potential_transfac_id = row['DBID.1']
            x = pd.read_csv(os.path.join(args.pwm_folder, pwm_file), sep="\t", index_col=0)
            if len(x) == 0:
                # Transfac2meme tool requires the $ character to be replaced with _
                transfac_id = potential_transfac_id.replace('$', '_')
                missing_motifs[transfac_id] = pwm_name
            else:
                motif_matrices[pwm_name] = np.array([row[1] for row in x.iterrows()])
                    
    except Exception as e:
        print e
        pass

print 'Made first pass through CisBP folder at {}'.format(args.pwm_folder)
print 'Found {} files, of which {} were empty (Assumed to be transfac)'.format(len(pwm_to_tf_names.keys()), len(missing_motifs.keys()))
print 'As a result, only {} PWM matrices could be made'.format(len(motif_matrices.keys()))
print 'TF distribution:'
print '{} motifs have 1 TF'.format(len([x for x in pwm_to_tf_names.keys() if len(pwm_to_tf_names[x]) == 1]))
print '{} motifs have 2 TFs'.format(len([x for x in pwm_to_tf_names.keys() if len(pwm_to_tf_names[x]) == 2]))
print '{} motifs have 3 TFs'.format(len([x for x in pwm_to_tf_names.keys() if len(pwm_to_tf_names[x]) == 3]))
print '{} motifs have 4 or more TFs'.format(len([x for x in pwm_to_tf_names.keys() if len(pwm_to_tf_names[x]) > 3]))
print ''

# Output the missing motifs
with open(os.path.join('transfac_ids'), 'w') as outfile:
    outfile.writelines('\n'.join(missing_motifs.keys()))

def convert_transfac_output(current_lines, current_motif):
    valid_lines = current_lines[2:-1] # lines in the pwm
    pwm_name = missing_motifs[current_motif]
    motif_matrices[pwm_name] = np.array([[float(y) for y in (x.lstrip(' ').rstrip('\n').split())] for x in valid_lines])

"""
# Iterate through transfac motifs, filling out those same data structures
with open(transfac_file) as f:
    current_motif = None
    current_lines = []
    try:
        for line in f.readlines():
            if line.startswith('MOTIF'):
                if current_motif:
                    # save the values from the PWM into the dictionary
                    convert_transfac_output(current_lines, current_motif)
                    current_lines = []
                current_motif = line.rstrip().lstrip('MOTIF').split(' ')[1]
                if args.debug:
                    print current_motif
            else:
                current_lines.append(line)
        # Handle very last motif in transfac file
        convert_transfac_output(current_lines, current_motif)

    except Exception as e:
        print e
"""

print 'After parsing the transfac file {}, the number of motif PWMs is {}'.format(transfac_file, len(motif_matrices.keys()))
print 'The following motifs were neither in transfac nor cisbp downloads: {}'.format(set(pwm_to_tf_names.keys())- set(motif_matrices.keys()))

if not args.mouse:

    cisbp_pwm_names = motif_matrices.keys()
    cisbp_tf_candidates = set([item for sublist in pwm_to_tf_names.values() for item in sublist])


    # Iterate through Encode motifs, filling out those same data structures
    encode_motifs = {}
    encode_extra_info = {}

    # Convert ENCODE motifs
    with open(encode_file) as f:
        current_motif = None
        current_lines = []
        try:
            for line in f.readlines():
                if line.startswith('>'):
                    if current_motif:
                        encode_motifs[current_motif] = current_lines
                        current_lines = []
                    current_motif = line.rstrip().lstrip('>')
                    if args.debug:
                        print current_motif
                else:
                    # Remove first character (A,C,T,G,R,Y, etc.)
                    current_lines.append(' '.join(line.split(' ')[1:]))
            encode_motifs[current_motif] = current_lines
        except Exception as e:
            print e
            pass

    # Loop through encode motifs to create numpy matrices, not string
    for key in encode_motifs:
        # parsing to deal with ENCODE's unique motif_info
        motif_info = key.split()
        pwm_name = motif_info[0]
        tf_name = pwm_name.split('_')[0]
        extra_info = motif_info[-1]
        motif_matrices[pwm_name] = np.array([[float(y) for y in x.rstrip('\n').split(' ')] for x in encode_motifs[key]])
        pwm_to_tf_names[pwm_name] = [tf_name]
        encode_extra_info[pwm_name] = extra_info
        
    encode_pwm_names = list(set(motif_matrices.keys()) - set(cisbp_pwm_names))

    print 'After parsing the encode file {}, the number of motif PWMs is {}'.format(encode_file, len(motif_matrices.keys()))
    print ''

    # Loop through encode motifs to find correlations against CisBP motifs
    # Or correlations against the inverse (updown flip)

    # np.corrcoef(motif_matrices[0].flatten(), np.flipud(motif_matrices[0]).flatten())

    num_correlations = 0
    threshold = .95

    duplicate_matrices = []
    duplicate_matrices_for_existing_tfs = []

    # Quick check: are all motifs unique? Defined as having a correlation lower than .95 (tho 2013 ENCODE used .75 as threshold)
    for i in cisbp_pwm_names:
        for j in encode_pwm_names:
            if len(motif_matrices[i]) == len(motif_matrices[j]):
                corrcoef = np.corrcoef(motif_matrices[i].flatten(), motif_matrices[j].flatten())
                if corrcoef[0][1] > threshold:
                    num_correlations = num_correlations + 1
                    # mark the encode motif for deletion
                    duplicate_matrices.append(j)
                    if pwm_to_tf_names[j][0] in pwm_to_tf_names[i]:
                        duplicate_matrices_for_existing_tfs.append(j)
                        if args.debug:
                            print 'Found correlation coefficient higher than {}, {}, at {} and {}'.format(threshold, corrcoef, i, j)
                            print 'TF1: {} vs TF2: {}'.format(pwm_to_tf_names[i], pwm_to_tf_names[j])

    print 'total correlations found in ENCODE vs 2016 dataset: {} with {} threshold'.format(num_correlations, threshold)
    print 'After parsing the encode file {}, the number of motif PWMs is {}'.format(encode_file, len(motif_matrices.keys()))
    print '{} encode motifs marked for deletion'.format(len(duplicate_matrices))
    print '{} encode motifs for existing TFs from CisBP marked for deletion'.format(len(set(duplicate_matrices_for_existing_tfs)))

    # delete dupes:
    for x in set(duplicate_matrices_for_existing_tfs):
        del motif_matrices[x]
        del pwm_to_tf_names[x]
    print 'After deleting encode duplicates the number of motif PWMs is {}'.format(len(motif_matrices.keys()))
    encode_pwm_names = list(set(motif_matrices.keys()) - set(cisbp_pwm_names))

    ## Second Section: assigning TFs to motifs

    all_gene_names = []
    with open(gene_list_file) as f:
        try:
            for line in f.readlines():
                attr = dict(item.strip().split(' ') for item in line.split('\t')[8].strip('\n').split(';') if item)
                all_gene_names.append(attr['gene_name'].replace('"', ''))
        except Exception as e:
            print e
            pass

    goa_tf_file = '/Users/ndeveaux/Data/Human_DNA_binding_proteins/human_tf_list_from_goa_2016_10_17.txt'
    goa_tfs = list(pd.read_csv(goa_tf_file, index_col=0).index)
    g = pd.read_csv("/Users/ndeveaux/Data/Human_DNA_binding_proteins/human_hgall_tf_list_from_tfclass_2016_10_17.txt")
    tf_class_tfs = [x.rstrip('\t') for x in list(g['HumanTF\tMouseTF'])]
    possible_human_tfs = set(tf_class_tfs).union(set(goa_tfs))
    existing_tf_candidates = set([item for sublist in [pwm_to_tf_names[k] for k in motif_matrices.keys()] for item in sublist])

    """
    potential_maps = {}
    for i in pwm_to_tf_names.values():
        for tf in i:
            if tf not in all_gene_names:
                potential_maps[tf] = []
                print 'Motif TF {} not in all gene names'.format(tf)
                for j in goa_tfs:
                    if j.startswith(tf):
                        potential_maps[tf].append(j)
                        print 'But it could be the prefix for {}'.format(j)

    # this is a flat value map
    potential_prefix_tf_candidates = set([item for sublist in potential_maps.values() for item in sublist])
    print 'Here are the TFs you could gain motifs for if you use prefix/generic motifs: {}'.format(potential_prefix_tf_candidates - existing_tf_candidates)"""

    covered_tfs = (existing_tf_candidates).intersection(possible_human_tfs)
    percentage = float(len(covered_tfs)) / float(len(possible_human_tfs))
    print 'TF coverage: {}, {} out of {} motif-tfs and {} annotated tfs and dna binding proteins'.format(percentage, len(covered_tfs), len(existing_tf_candidates), len(possible_human_tfs))
    really_all_gene_names = set(all_gene_names).union(possible_human_tfs)
    print 'TF coverage: {} possible tfs from all {} human genes'.format(len(really_all_gene_names.intersection(existing_tf_candidates)) , len(really_all_gene_names))

    # Third Section: write to output:

    for key in encode_pwm_names:
        # import pdb; pdb.set_trace()
        df2 = pd.DataFrame([[key, pwm_to_tf_names[key][0], encode_extra_info[key]]], columns=['Motif_ID', 'TF_Name', 'MSource_Type'])
        tf_information = tf_information.append(df2)
    
    # Replace all NaN's added from the encode data appending
    tf_information = tf_information.replace(np.nan,'', regex=True)
    
    # Copy the TFInformation.txt input into a new TFInformation.txt:
    tf_information.to_csv("new_TF_information.txt", sep="\t", index=False)

    # Additional output file
    tf_to_motifs = {}
    for motif_name in motif_matrices.keys():
        tfs = pwm_to_tf_names[motif_name]
        for tf in tfs:
            if tf not in tf_to_motifs.keys():
                tf_to_motifs[tf] = [motif_name]
            else:
                tf_to_motifs[tf].append(motif_name)

    # print list of known human tfs, and motifs
    with open('human_tfs_to_motifs.txt', 'w') as outfile:
        lines = '\n'. join(['{}\t{}'.format(key, ','.join(val)) for key, val in tf_to_motifs.items() if key in possible_human_tfs])
        outfile.writelines(lines)

else:
    mouse_tf_file_curated_by_Emily = '/Users/ndeveaux/Data/Mouse_DNA_binding_proteins/DNA_binding_final4Matlab.txt'
    import pdb; pdb.set_trace()
    possible_mouse_tfs = list(pd.read_csv(mouse_tf_file_curated_by_Emily, index_col=0).index)

        # Additional output file
    tf_to_motifs = {}
    for motif_name in motif_matrices.keys():
        tfs = pwm_to_tf_names[motif_name]
        for tf in tfs:
            if tf not in tf_to_motifs.keys():
                tf_to_motifs[tf] = [motif_name]
            else:
                tf_to_motifs[tf].append(motif_name)

    # print list of known human tfs, and motifs
    with open('mouse_tfs_to_motifs.txt', 'w') as outfile:
        lines = '\n'. join(['{}\t{}'.format(key, ','.join(val)) for key, val in tf_to_motifs.items() if key in possible_mouse_tfs])
        outfile.writelines(lines)




for motif in motif_matrices:
    matrix = motif_matrices[motif]
    tf_names = ','.join(pwm_to_tf_names[motif])
    result = template.substitute(name=motif,
                                            tf_name=tf_names,
                                             width=len(matrix),
                                             nsites=len(matrix))
    result += "\n".join(["  ".join([ "{:.6}".format(x) for x in y]) for y in matrix]) + '\n'
    output += result 

# Write output meme file:
with open(os.path.join(args.out_file), 'w') as outfile:
    outfile.writelines(output)

