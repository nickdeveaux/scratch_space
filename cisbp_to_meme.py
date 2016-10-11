#!/usr/bin/python

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

args = parser.parse_args()

with open(os.path.join(pwm.pwm_dir(), "cisbp_header.txt")) as input:
    header = input.read()

output = header

with open(os.path.join(pwm.pwm_dir(), "cisbp_template.txt")) as input:
    template = MyTemplate(input.read())

if args.information_file:
    tf_information = pd.read_csv(args.information_file, sep="\t", header=0)

new_tf_names = {}
cisbp_motif_matrices = {}

# Due to licensing issues, the TRANSFAC motifs (that we have a license for!) are not able to be downloaded from CisBP
# The solution is to save their IDs to a file that will be based in via the -id flag to Meme Suite's transfac2meme tool
missing_motifs = {}

for pwm_file in os.listdir(args.pwm_folder):
    try:
        pwm_name = ".".join(pwm_file.split(".")[:-1])
        if pwm_name in tf_information['Motif_ID'].tolist():
            print pwm_name
            index = tf_information['Motif_ID'].tolist().index(pwm_name)
            tf_name = tf_information.loc[index,'TF_Name']
            new_tf_names[pwm_name] = tf_name
            ens_id = tf_information.loc[index,'DBID']
            x = pd.read_csv(os.path.join(args.pwm_folder, pwm_file), sep="\t", index_col=0)
            if len(x) == 0:
                transfac_id = tf_information.loc[index,'DBID.1']
                # Transfac2meme tool requires the $ character to be replaced with _
                transfac_id = transfac_id.replace('$', '_')
                missing_motifs[transfac_id] = [pwm_name, tf_name, ens_id]
            else:
                result = template.substitute(name=pwm_name,
                                            tf_name=tf_name,
                                            ens_id=ens_id,
                                             width=len(x),
                                             nsites=len(x))
                    
                for row in x.iterrows():
                    result += "  ".join(row[1].map(lambda x: "{:.6}".format(x))) + "\n"
    
                cisbp_motif_matrices[pwm_name] = np.array([row[1] for row in x.iterrows()])
        
                output += result
            
    except Exception as e:
        print e
        pass

transfac_motifs = {}

# Parse the transfac output:
with open("/Users/ndeveaux/Data/motifs/transfac_output.meme") as f:
    current_motif = None
    current_lines = []
    try:
        for line in f.readlines():
            if line.startswith('MOTIF'):
                if current_motif:
                    # save the values from the PWM into the dictionary
                    valid_lines = current_lines[2:-1] # lines in the pwm
                    transfac_motifs[current_motif] = ''.join(valid_lines).rstrip('\n')
                    pwm_name, tf_name, ens_id = missing_motifs[current_motif]
                    result = template.substitute(name=pwm_name,
                                            tf_name=tf_name,
                                            ens_id=ens_id,
                                             width=len(valid_lines),
                                             nsites=len(valid_lines))
                    result += transfac_motifs[current_motif] + '\n'
                    output += result 
                    current_lines = []
                current_motif = line.rstrip().lstrip('MOTIF').split(' ')[1]
                print current_motif
            else:
                current_lines.append(line)
    except Exception as e:
        print e

# parse transfac into 
for i in transfac_motifs.keys():
    pwm_name, tf_name, ens_id = missing_motifs[i]
    cisbp_motif_matrices[pwm_name] = np.array([[float(y) for y in (x.lstrip(' ').split())] for x in transfac_motifs[i].split('\n')])


encode_motifs = {}

# Convert ENCODE motifs
with open("/Users/ndeveaux/Data/motifs/Encode_Site/motifs.txt") as f:
    current_motif = None
    current_lines = []
    try:
        for line in f.readlines():
            if line.startswith('>'):
                if current_motif:
                    encode_motifs[current_motif] = current_lines
                    current_lines = []
                current_motif = line.rstrip().lstrip('>')
                print current_motif
            else:
                # Remove first character (A,C,T,G,R,Y, etc.)
                current_lines.append(' '.join(line.split(' ')[1:]))
    except Exception as e:
        print e
        pass

encode_motif_matrices = []
for val in encode_motifs.values():
    encode_motif_matrices.append(np.array([[float(y) for y in x.rstrip('\n').split(' ')] for x in val]))
    
# Loop through encode motifs to create numpy matrices, not string


# Loop through encode motifs to find correlations against CisBP motifs
# Or correlations against the inverse (updown flip)

# np.corrcoef(motif_matrices[0].flatten(), np.flipud(motif_matrices[0]).flatten())

num_correlations = 0
threshold = .95
key_list = list(encode_motifs.keys())
# Quick check: are all motifs unique? Defined as having a correlation lower than .98 (tho 2013 ENCODE used .75 as threshold)
for i in cisbp_motif_matrices.keys():
    for j in range(len(encode_motif_matrices)):
        if len(cisbp_motif_matrices[i]) == len(encode_motif_matrices[j]):
            corrcoef = np.corrcoef(cisbp_motif_matrices[i].flatten(), encode_motif_matrices[j].flatten())
            if corrcoef[0][1] > threshold:
                print 'Found correlation coefficient higher than {}, {}, at {} and {}'.format(threshold, corrcoef, i, key_list[j])
                print 'TF1: {} vs TF2: {}'.format(new_tf_names[i], key_list[j])
                num_correlations = num_correlations + 1
                if key_list[j] in encode_motifs.keys():
                    del encode_motifs[key_list[j]]

print 'total correlations found in ENCODE vs 2016 dataset: {} with {} threshold'.format(num_correlations, threshold)

print len(encode_motifs)
print len(encode_motif_matrices)

# loop through remaining encode motifs and turn them into text via template:
output += '\n'
for i in encode_motifs:
    try:
        motif_info = i.split() # split without any arguments handles tabs and spaces simultaneously
        pwm_name = motif_info[0]
        extra_info = motif_info[-1]
        tf_name = pwm_name.split('_')[0]
        result = template.substitute(name=pwm_name,
                                            tf_name=tf_name,
                                            ens_id=extra_info,
                                             width=len(encode_motifs[i]),
                                             nsites=len(encode_motifs[i]))
        for j in encode_motifs[i]:
            result += "  ".join([ "{:.6}".format(x) for x in j.rstrip('\n').split(' ')]) + "\n"
        output += result
    except Exception as e:
        print e
        print i
        pass

import pdb; pdb.set_trace()
with open(os.path.join(args.out_file), 'w') as outfile:
    outfile.writelines(output)

with open(os.path.join('transfac_ids'), 'w') as outfile:
    outfile.writelines('\n'.join(missing_motifs.keys()))


"""
missing_tfs = []
missing_motifs = []
composite_motif_lines = []
is_composite_motif = False
all_old_motifs = []
all_old_tfs = []

# Add existing data from hg19_em.meme
with open("/Users/ndeveaux/Data/motifs/hg19_em.meme") as f:
    for line in f.readlines():
        if line.startswith('MOTIF'):
            if len(line.rstrip().split(' ')) < 3:
                print line
            else:
                is_composite_motif = False
                motif, tf = line.rstrip().split(' ')[1:3]
                all_old_motifs.append(motif.split('.')[0])
                all_old_tfs.append(tf)
                if '::' in tf:
                    is_composite_motif = True
                if not tf in new_tf_names:
                    missing_tfs.append(tf)
                if not motif.split('.')[0] in new_motif_names:
                    missing_motifs.append(motif.split('.')[0])
        if is_composite_motif:
            composite_motif_lines.append(line)


print 'missing motifs: {}, ex: {}'.format(str(len(set(missing_motifs))), str(missing_motifs[1:5]))
print 'missing tfs: {}, ex: {}'.format(str(len(set(missing_tfs))), str(missing_tfs[1:5]))
composite_subset = [x for x in missing_tfs if '::' in x]
print 'missing tfs with composites: {}, ex: {}'.format(str(len(set(composite_subset))), str(composite_subset[1:5]))

new_motifs = set(new_motif_names) - set(all_old_motifs)
new_tfs = set(new_tf_names) - set(all_old_tfs)

print 'new motifs: {}, ex:{}'.format(str(len(new_motifs)), str(list(new_motifs)[1:5]))
print 'new tfs: {}, ex:{}'.format(str(len(new_tfs)), str(list(new_tfs)[1:5]))

print 'extra motifs in TF_Information.txt: {}'.format(str(len(set([x for x in tf_information['Motif_ID'].tolist() if x not in new_motif_names]))))
print 'extra tfs in TF_Information.txt: {}'.format(str(len(set([x for x in tf_information['TF_Name'].tolist() if x not in new_tf_names]))))

print 'total size, in motif count, of old hg19_em.meme motif db: {}'.format(len(all_old_motifs))
print 'total size, in motif count, of new {} motif db: {}'.format(args.information_file, len(new_motif_names))
"""
