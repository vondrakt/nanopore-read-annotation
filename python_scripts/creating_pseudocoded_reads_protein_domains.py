#! /usr/bin/env python3

from optparse import OptionParser
import re
from Bio import SeqIO

# THE READS ARE DELIVERED AS A STANDARD INPUT
# AND THE OUTPUT IS DELIVERED AS A STANDARD OUTPUT


# defining the arguments that need to be passed to the script
arguments = OptionParser()

arguments.add_option('-i', '--reads', dest='reads', help='pseudocoded_reads')
arguments.add_option('-c', '--coding', dest='coding_table', help='input the coding table')
arguments.add_option('-g', '--gff', dest='gff_file', help='input the gff file')

(options, args) = arguments.parse_args()
if options.coding_table is None or options.gff_file is None or options.reads is None:
    # if one of the arguments is missing
    print('\n----------> A mandatory option is missing !\n')  # raise an error
    arguments.print_help()  # and provide the help
    exit(-1)  # exit the script

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script takes the existing pseudocoded reads, created based on the lastz tabular output and adds annotation for 
# protein domains into the pseudocoded reads.
# The annotation for the domains is kept in the gff file.
# In order to match the gff hits to the corresponding pseudocode, a coding table is needed.

dictionary_upper = {}  # contains the pseudocodes of + oriented repeats
dictionary_lower = {}  # contains pseudocodes of - oriented repeats

# parsing the coding table
# opening the file with the coding table
with open(options.coding_table) as f: 
    # iterating over lines in file
    for line in f:
        # the lines are split into individual items
        items = line.split()

        # if the pseudocode is defined as forward in the coding table 
        if items[1] == 'F': 
            # store it in the dictionary_upper
            dictionary_upper[items[0]] = items[2]
        # if not
        else:
            # store it in the dictionary_lower variable
            dictionary_lower[items[0]] = items[2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dictionary_reads = {}  # contains the reads submitted for the analysis

# opening the file with the reads
with open(options.reads) as f:

    # iterating over all reads
    for read in SeqIO.parse(f, 'fasta'):

        # storing it within the dictionary_reads
        dictionary_reads[read.id] = str(read.seq) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# opening the gff file that contains the coordinates of protein domains
with open(options.gff_file) as f:
    # iterating over lines of the gff file
    for line in f:
        # if the line is empty
        if line == '\n':
            # continue to the next line
            continue

        # if the line begins with a '#'
        if line[0] == '#':
            # continue to the next line
            continue

        # split each line into separate items based on the tab
        items = line.split()

        # if the read ID from the gff exists in the dictionary_reads variable
        if items[0] in dictionary_reads.keys():

            # find the classification of the protein domain
            hit = re.search('Final_Classification=[\w\|\-/]+', items[8]).group()
            classification = re.sub('Final_Classification=', '', hit)
            hit = re.search('Name=\w+', items[8]).group()
            domain = re.sub('Name=', '', hit)

            # the variable code, stores the final classification of the domain
            code = classification + '__' + domain


            # deciding on orientation, if the orientation is '-'
            if items[6] == '-':
                
                # check if the lineage of the protein domain is recorded in the coding table
                if code in dictionary_lower.keys():
                    # if yes store the pseudocode in the variable pseudocode
                    pseudocode = dictionary_lower[code]

                # if it is not
                else:
                    # the pseudocode to be added to the read is 'z'
                    pseudocode = 'z'

            # if the orientation is '+'
            if items[6] == '+':
                # check if the lineage of the protein domain is recorded in the coding table
                if code in dictionary_upper.keys():
                    pseudocode = dictionary_upper[code]

                # if it is not
                else:
                    # the pseudocode to be added to the read is 'Z'
                    pseudocode = 'Z'

            # changing that part of the read
            # calculating the position to which the pseudocode string will be added to the read
            start = int(items[3]) - 1
            end = int(items[4]) - 1

            # turning the read string into a list
            read_list = list(dictionary_reads[items[0]])
            # variable i will be used as a counter to the range where the read will be changed
            i = start

            # while i is smaller than or equal to the ending position
            while i <= end:

                # if the variable i is smaller than the last position of the read
                if i > len(read_list)-1:
                    print(items[0])

		# if the position in the read is equal to '0'
                if read_list[i] == '0':
                    # change that position to the pseudocode
                    read_list[i] = pseudocode
                # increase the i variable
                i += 1

            # join the read list into a string and reassign it to the corresponding key in the dictionary_reads variable
            dictionary_reads[items[0]] = ''.join(read_list)

# for the items in the dictionary_reads variable
for key, val in dictionary_reads.items():
    # print the key and value in the fasta format
    print('>'+key+'\n'+val)
