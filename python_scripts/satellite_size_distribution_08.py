#!/usr/bin/env python3

from optparse import OptionParser
from Bio import SeqIO

# defining the arguments which can be passed to the script
arguments = OptionParser()

arguments.add_option('-i', '--in', dest='reads', help='text file with pseudocoded reads')
arguments.add_option('-s', '--start', dest='start_pos', help='minimum position from which read are truncated')
arguments.add_option('-o', '--out', dest='output', help='file to which to write output')
arguments.add_option('-c', '--coding', dest='coding_table', help='file of the coding table')
(options, args) = arguments.parse_args()
if (options.reads is None or options.output is None or options.start_pos is None or options.coding_table is None):  # if one of the arguments is not provided
        print('\n----------> A mandatory option is missing !\n')   # raise an error

        arguments.print_help()  # and provide the help
        exit(-1)  # exit the script

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the coding_table dictionary will have the keys as the main satellite name
# and the value will be the pseudocode assigned to the satellite name
coding_table = {}

# open the coding table file
with open(options.coding_table) as c:
    # iterate over all lines
    for line in c:
        # split columns of each line in separate items
        items = line.split()
        # assign the pseudocodes as keys of the dictionary and values as the satellite names
        coding_table[items[2]] = items[0]

# because X and 0 are not defined in the coding table but do exist in pseudocoded reads
# these two entries need to be added to the coding_table variable
coding_table['X'] = 'X'
coding_table['0'] = '0'

# the dictionary dictionary_reads will store the read IDs as keys and pseudocoded strings as values
dictionary_reads = {}

# opening the read file
with open(options.reads) as f:
    # iterating over sequences in multifasta file
    for seq in SeqIO.parse(f, 'fasta'):
        # adding entries to the dictionary
        dictionary_reads[seq.id] = str(seq.seq)

# length_data is going to be a list of lists which will store information about
# individual arrays of satellites, their length, whether they are intact or truncated anthe read ID
# and read length where the arrays was found
length_data = []

# iterating over items of dictionary
for key, val in dictionary_reads.items():

    # indexing each string position of the pseudocode string
    for idx, n in enumerate(val):

        # check if the list is empty and if it is, add the first value
        if len(length_data) == 0:
                # if there is a part of the pseudocode that is not in the coding table
                if n not in coding_table.keys():
                    # print that code and return the message
                    print(n)
                    print('code not in coding table')
                # get the satellite name of the pseudocode array from the coding_table
                sat_name = coding_table[n]
                # create a list with the read ID, first and second index of the array, satellite name and pseudocode
                length_row = [key, idx, idx, sat_name, n]
                # append to existing list of lists
                length_data.append(length_row)

        # if the length_data is not empty and the position has the same pseudocode as the last pseudocode within the length_data
        if length_data[-1][4] == n:

                # check if the array is continuous
                if idx - length_data[-1][2] == 1 or idx - length_data[-1][2] == 0:
                    # if it is then increment the length of the bed by one position
                    length_data[-1][2] = idx

                # if the bed is not continuous
                else:
                    # check if the pseudocode is not in the coding table
                    if n not in coding_table.keys():
                        print(n)
                        print('code not in coding table')

                    # get the satellite name of the pseudocode array from the coding_table
                    sat_name = coding_table[n]
                    # create a new entry for the length_data with a new interval
                    length_row = [key, idx, idx, sat_name, n]
                    # add to length_data
                    length_data.append(length_row)

        # if the pseudocodes are different then start new interval and repeat the steps described above
        # for adding a new entry
        else:
                if n not in coding_table.keys():
                    print(n)
                    print('code not in coding table')
                sat_name = coding_table[n]
                length_row = [key, idx, idx, sat_name, n]
                length_data.append(length_row)

# opening the output file
out = open(options.output, 'w')
# iterating over length_data
for i in length_data:
    # store the read length
    read_length = len(dictionary_reads[i[0]])
    # calculate the length of the array
    length = i[2] - i[1] + 1
    # check if the array start below the first limited number of bases or after them
    if i[1] < int(options.start_pos) or i[2] > (read_length - int(options.start_pos)):
        # if it does then the array is truncated
        # create the correct output format and write out line
        line = str(i[3])+'\t'+str(length)+'\t'+str(i[4])+'\t'+str(read_length)+'\t'+'T'+'\t'+i[0]+'\n'
    else:
        # if it does not then the array is intact
        # create the correct output format and write out line
        line = str(i[3])+'\t'+str(length)+'\t'+str(i[4])+'\t'+str(read_length)+'\t'+'I'+'\t'+i[0]+'\n'

    out.write(line)






