#!/usr/bin/env python3

from optparse import OptionParser
import re

# modified version of the script where the ends of the reads are cut off due to lower quality

# defining the arguments which can be passed to the script
arguments = OptionParser()

arguments.add_option('-r', '--reads', dest='reads', help='text file with pseudocoded reads')
arguments.add_option('-w', '--window', dest='window', help='size of window')
arguments.add_option('-c', '--coding', dest='coding', help='coding table')
arguments.add_option('-s', '--start', dest='start', help='starting position')
arguments.add_option('-o', '--out',  dest='output', help='name of output file')

(options, args) = arguments.parse_args()
if options.reads is None or options.window is None or options.coding is None\
        or options.output is None or options.start is None:  # if one of the arguments is not provided
        print('\n----------> A mandatory option is missing !\n')   # raise an error

        arguments.print_help()  # and provide the help
        exit(-1)  # exit the script

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# creating the dictionary which stores the pseudocodes and their complements
complement_dictionary = {}
# the list forward will store all the pseudocodes for the forward orientation
forward = []
# the list reverse will store all the pseudocodes for the reverse orientation
reverse = []

# opening the coding table file
with open(options.coding) as c:
    # iterating over each line in file
    for line in c:
        # dividing columns of each line in file
        items = line.split()
        # if the pseudocode is classified as forward
        if items[1] == 'F':
            # append it to the forward list
            forward.append(items[2])
        # if the pseudocode is classified as reverse
        else:
            # append it to the reverse list
            reverse.append(items[2])

# variable i will be used for counting the length of the lists
i = 0
# while i is smaller than the length of the list forward
while i < len(forward):
    # add the complements to the complement_dictionary
    complement_dictionary[forward[i]] = reverse[i]
    complement_dictionary[reverse[i]] = forward[i]
    # increment i
    i += 1

# because X and 0 are not defined in the coding table but do exist in pseudocoded reads
# these two entries need to be added to the complement_dictionary
complement_dictionary['0'] = '0'
complement_dictionary['X'] = 'X'
complement_dictionary['Z'] = 'z'
complement_dictionary['z'] = 'Z'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to create the reverse complement of the read
def complement(read, dictionary):
    # the complement_read list will store the reverse complement of the read being processed
    complement_read = []
    # the counter will start from the last position
    i = -1

    # while i is larger than the negative length of the read
    while i > -len(read) - 1:
        # the complement of the current position is appended to the list
        complement_read.append(dictionary[read[i]])
        # decrease i by -1
        i -= 1
    # join the complement read list into a string
    complement_read = ''.join(complement_read)

    # return the complement_read
    return complement_read
    # the output of the function is a string


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function which creates the left and right dictionaries
# the profiles of the neighborhood will be kept as dictionaries with keys for each satellite group
# values of these dictionaries will be lists of 0 the length of the window given in the options

def create_left_and_right(coding_table, window_size):
    # creating left and right dictionaries
    left = {}
    right = {}

    # opening the file containing the coding table
    with open(coding_table) as c:
        # reading lines of file
        for line in c:
            # splitting into a list
            items = line.split()
            # adding entries to both dictionaries with keys from the coding table and lists of 0
            left[items[2]] = list(
                map(int, '0' * int(window_size)))
            right[items[2]] = list(map(int, '0' * int(window_size)))

    # adding entries that are not within the coding table
    left['0'] = list(map(int, '0' * int(window_size)))
    right['0'] = list(map(int, '0' * int(window_size)))

    left['X'] = list(map(int, '0' * int(window_size)))
    right['X'] = list(map(int, '0' * int(window_size)))

    left['Z'] = list(map(int, '0' * int(window_size)))
    right['Z'] = list(map(int, '0' * int(window_size)))

    left['z'] = list(map(int, '0' * int(window_size)))
    right['z'] = list(map(int, '0' * int(window_size)))

    # return the two dictionaries
    return left, right



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function which will find all the left neighbours of a pseudocode array and increment the dictionary
def left_neighbour(pseudo, working_read, left):
    # variable first stores the first occurence of the pseudocode within the read
    first = None
    # list first_positions stores all the first positions within the read of the pseudocode being processed
    first_positions = []
    # iterating over all values and indices
    for idx, val in enumerate(working_read):
        # if the value of string is different from the pseudocode
        if val != pseudo:
            # first is None
            first = None
        # if the value of string is the same as the pseudocode
        else:
            # check if it is the first occurence
            if first is None:
                # if it is set the first to the index of the first occurence
                first = idx
                # append it to the first_positions list
                first_positions.append(first)
            # if it is not the first occurence continue to next
            else:
                continue

    # extracting the windows sequences
    for i in first_positions:

        # calculating the starting position of the window
        window_start = i - int(options.window)

        # if it starts below the boundaries of the read
        if window_start < 0:
            # start from 0
            window_seq = working_read[:i]
        # if it does not
        else:
            # start from the given position
            window_seq = working_read[window_start:i]

        # iterating over the window sequence
        i = -1

        # while i is larger than the negative length of the window
        while i > -len(window_seq) - 1:
            # increment the window by 1
            left[window_seq[i]][i] += 1
            # decrease i by 1
            i -= 1

    # return the left dictionary
    return left

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function which will find all ending positions of a pseudocode array and increment the dictionary
def right_neighbour(pseudo, working_read, right):
    # list last_positions stores all the last positions within the read of the pseudocode being processed
    last_positions = []
    # variable last stores the last occurence of the pseudocode within the read
    last = None
    i = 0

    # while i is smaller than the length of the pseudocoded read
    while i < len(working_read) - 1:

        # check if the position is different to the pseudocode
        if working_read[i] != pseudo:
            # if it is last is None
            last = None
            # increment i by 1
            i += 1

        # if the positions is the same as the pseudocode
        else:

            # while the positions are the same as the pseudocode
            while working_read[i] == pseudo:

                # increment i by 1
                i += 1

                # if i is the same as the length of the read
                if i == len(working_read) - 1:
                    # break the loop
                    break

            # set last to the last i position
            last = i
            # append the position to last_positions
            last_positions.append(last)

    # calculating the end positions of the window
    # iterating over last positions
    for i in last_positions:
        # calculate the ending of the window
        window_end = i + int(options.window)

        # if the end is beyond the length of the read
        if window_end > len(working_read)-1:
            # stop the window sequence at the read end
            window_seq = working_read[i:]

            # if the length of the window sequence is 1 or 0
            if len(window_seq) <= 1:
                # conitnue to next window
                continue

        # if the end is within the read length
        else:
            # extract the window sequence
            window_seq = working_read[i:window_end]

        # iterate over positions within the window string
        for idx, val in enumerate(window_seq):
            # increment the positions in the dictionary by 1
            right[val][idx] += 1

    # return the right dictionary
    return right


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function which creates the structure of the output table and writes it to file
def write_table(window_size, right, left, output_name, pseudo):
    # the first line of the file will be the positions from 0 to the window length
    # both on the right and left side
    min_position = -int(window_size)
    max_position = int(window_size)
    # variable with the range of bases which will be the first line in the file
    base_range = list(range(min_position, max_position + 1))
    # converting base_range list to a tab separated string
    base_range = '\t'.join(str(base_range).split(','))
    base_range = re.sub('\[', '', base_range)
    base_range = re.sub('\]', '', base_range)

    # opening the output file
    # the output file will have the name of the pseudocode in the file name
    out_name = output_name + '_' + pseudo
    out = open(out_name, 'w')

    # first write the base range to the file
    out.write('\t' + base_range + '\n')

    for key in right.keys():
        # this line just coverts the numeric list to a tab separated string
        row_of_file_left = '\t'.join(str(left[key]).split(','))
        row_of_file_left = re.sub('\[', '', row_of_file_left)
        row_of_file_left = re.sub('\]', '', row_of_file_left)
        row_of_file_right = '\t'.join(str(right[key]).split(','))
        row_of_file_right = re.sub('\[', '', row_of_file_right)
        row_of_file_right = re.sub('\]', '', row_of_file_right)

        # creating the line format
        line = key + '\t' + row_of_file_left + '\t' + '0' + '\t' + row_of_file_right
        # now we write the lines to the output file
        out.write(line + '\n')
    # close the output file
    out.close()
    # message is printed every time the output for a pseudocode group is written
    print('Output table has been written')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function which will sum up all the bases that have been counted
def count_bases(right, left, window, pseudo, output_name):
    # list for the sum of the bases on the right
    bases_on_right = []
    # list for the sum of the bases on the left
    bases_on_left = []

    i = 0
    # while i is not the length of the window
    while i != int(window):
        # set the r and l variables to 0
        # they will count the same by column in the dictionaries
        r = 0
        l = 0

        # iterating over all keys on the right
        for key in right.keys():
            # iterate through all the keys for just the position i
            # increment the r and l variables
            r = r + right[key][i]
            l = l + left[key][i]

        # append the sum of the position to the corresponding lists
        bases_on_right.append(r)
        bases_on_left.append(l)
        # increment i by 1
        i += 1

    # the 0 will represent the collapse satellite arrays and will divide the left and right side
    bases_on_left.append(0)

    # all_bases is a list that combines the sums of the left and right sides
    all_bases = bases_on_left + bases_on_right

    # the first line of the file will be the positions
    min_position = -int(window)
    max_position = int(window)
    # variable with the range of bases which will be the first line in the file
    base_range = list(range(min_position, max_position + 1))
    # converting base_range list to a tab separated string
    base_range = '\t'.join(str(base_range).split(','))
    base_range = re.sub('\[', '', base_range)
    base_range = re.sub('\]', '', base_range)

    # opening the output file
    # the output file will have the name of the pseudocode in the file name
    out_name = 'base_count_' + pseudo + '_' + output_name
    out = open(out_name, 'w')

    # join the list into a string divided by \t
    all_bases = '\t'.join(str(all_bases).split(','))
    all_bases = re.sub('\[', '', all_bases)
    all_bases = re.sub('\]', '', all_bases)

    # first write the base range to the file
    out.write(base_range + '\n')
    out.write(all_bases)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# getting the list of pseudocodes that need to be processed
pseudocodes = []

# open the coding table file to get the pseudocodes
with open(options.coding) as c:
    # iterate over all lines
    for line in c:
        # divide the columns of each line into separate
        items = line.split()

        # only append those pseudocodes which are classified as forward
        if items[1] == 'F':
            pseudocodes.append(items[2])

# append pseudocode Z which represents all non-Ogre domains found in reads
pseudocodes.append('Z')
pseudocodes_set = set(pseudocodes)
pseudocodes = list(pseudocodes_set)

# iterate over pseudocodes
for i in pseudocodes:
    # for each individual pseudocode create a left and right dictionary
    left, right = create_left_and_right(options.coding, options.window)

    # open the reads file
    with open(options.reads) as r:
        # iterate over all lines
        for line in r:

            # if the line is empty, continue to next one
            if line == '\n':
                continue

            # if the line starts with >, conitnue to next one
            if line[0] == '>':
                continue

            # first analyse the forward oriented read
            # remove the \n
            read = re.sub('\n', '', line)
            # extract the part of the reads that falls in between the first number of bases given in the options
            end = len(read) - int(options.start) - 1
            read = read[int(options.start):end]
            # increment the left neighbour dictionary
            left = left_neighbour(i, read, left)
            # increment the right neighbour dictionary
            right = right_neighbour(i, read, right)

            # now analyse the reverse complement
            reverse_complement = complement(read, complement_dictionary)
            # extract the part of the reads that falls in between the first number of bases given in the options
            end = len(reverse_complement) - int(options.start) - 1
            reverse_complement = reverse_complement[int(options.start):end]
            # increment the left neighbour dictionary
            left = left_neighbour(i, reverse_complement, left)
            # increment the right neighbour dictionary
            right = right_neighbour(i, reverse_complement, right)

    # write the incrementation to the output
    write_table(options.window, right, left, options.output, i)
    # count the bases processed and write to output
    count_bases(right, left, options.window, i, options.output)
