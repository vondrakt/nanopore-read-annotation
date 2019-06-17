#! /usr/bin/env python3


# HERE I CHECK THE PRIORITIES AS WELL AS IF THE ELEMENTS HAVE THE SAME PRIORITIES AND IF THEY ARE THE SAME NAME!


from optparse import OptionParser
import re
import sys

# defining the arguments that need to be passed to the script
arguments = OptionParser()


arguments.add_option('-c', '--coding', dest='coding_table', help='input the coding table')

(options, args) = arguments.parse_args()
if options.coding_table is None:
    # if one of the arguments is missing
    print('\n----------> A mandatory option is missing !\n')  # raise an error
    arguments.print_help()  # and provide the help
    exit(-1)  # exit the script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dictionary_upper = {}  # contains the pseudocodes of + oriented repeats
dictionary_lower = {}  # contains pseudocodes of - oriented repeats
dictionary_lengths = {}  # contains minimum allowed lengths for each repeat
dictionary_priorities = {}  # contains levels of priority for each element

# parsing the coding table

with open(options.coding_table) as f:  # opening the file with the coding table
    for line in f:  # iterating over lines in file
        items = line.split('\t')  # the lines are split into individual items
        items[4] = re.sub('\n', '', items[4])
        if items[1] == 'F':  # if the repeat is + oriented
            dictionary_upper[items[0]] = items[2]  # save the upper pseudocode
            dictionary_lengths[items[2]] = items[3]  # the minimal length for upper pseudocode
            dictionary_priorities[items[2]] = items[4]  # the priority for upper pseudocode
        else:  # if the repeat is - oriented
            dictionary_lower[items[0]] = items[2]  # save the lower pseudocode
            dictionary_lengths[items[2]] = items[3] # for lower pseudocode
            dictionary_priorities[items[2]] = items[4]  # yje priority for lower pseudocode


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# defining a dictionary which will have pseudocodes as keys and strings of 0 as values
# the strings of this dictionary will be incremented by 1 if a hit is present

def new_dictionary(upper_coding_table, lower_coding_table, read_length):
    d = {}  # dictionary which will contain the strings
    for key in upper_coding_table.keys():  # iterating over both dictionaries because the
        # repeat names as values are the same in both the dictionaries
        d[upper_coding_table[key]] = list(map(int, '0'*int(read_length)))
        d[lower_coding_table[key]] = list(map(int, '0'*int(read_length)))
    return d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# first the dictionary will be incremented based on the positions of the hit

def increment(dict, start, length, orientation, repeat_name, upper_coding_table, lower_coding_table):
    start = int(start) - 1  # calculating start
    end = start + int(length)  # calculating end

    if orientation == '+':  # if the hit is + oriented
        pseudocode = upper_coding_table[repeat_name]  # extracting the key of dict from upper coding table
        dict[pseudocode][start:end] = [i + 1 for i in dict[pseudocode][start:end]]   # incrementing based on hit

    else:  # if the hit is - oriented
        pseudocode = lower_coding_table[repeat_name]  # extracting the key of dict from lower coding table
        dict[pseudocode][start:end] = [i + 1 for i in dict[pseudocode][start:end]]  # incrementing based on hit

    return dict  # return the incremented dictionary

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# now the arrays which are shorter than the minimum length required by the coding table
# are filtered out

def filter(dict, length_table):
    for key, val in dict.items(): # iterating over keys and values of dictionary
        le = 0  # variable which will store the length of arrays
        min_length = int(length_table[key])  # variable which stores the minimum length required
        # for the arrays to be kept
        index = -1  # variable which stores the index of the string
        first = None  # variable which stores the first position of an array

        for n in val:  # iterating over values within string

            index += 1  # increment the index value for each position
            if n != 0:  # if the position is different from 0
                le += 1  # increment length by 1
                if not first:  # and if the first position is not defined
                    first = index  # define it as the first index where the position is different from 0

            if n == 0:  # if the position is 0

                if le < min_length and first:  # check if the array length is smaller than the minimum length
                    for ii in range(first, index):  # change this interval to 0
                        val[ii] = 0

                    le = 0  # set the length to 0 for the next array
                    first = None  # set the first position to None for next array

                if le > min_length:  #if the length is larger than minimum length
                    le = 0  # just set these variables to 0 for the next array
                    first = None

        if index == len(val) - 1:  # if the loop is at the end of the string
            if le < min_length and first:  # check the length of the last array
                for ii in range(first, index+1):
                    val[ii] = 0
    return dict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# after filtering short arrays, the dictionary must be collapsed into one pseudocoded sequences which
# will be the output of the script
# while creating the final pseudocode the priorities must be taken into consideration

def pseudocode_string(dict, priority_table, read_length):
    string = list(map(int, '0'*int(read_length)))  # assigning a string of 0 length of the read
    i = 0  # variable which stores the position within the string

    while i != read_length:  # while loop which iterates over positions in dictionary d

        position_pseudocode = 0   # for the first position the position_pseudocode and position_priority
        # are not defined
        position_priority = 0

        for key in dict:  # iterating over keys within dictionary

            if dict[key][i] != 0:

                if int(priority_table[key]) == position_priority:
                    position_pseudocode = 'X'
                    position_priority = int(priority_table[key])

                if int(priority_table[key]) > position_priority:
                    position_priority = int(priority_table[key])
                    position_pseudocode = key


        # the string is changed to match the pseudocode which had the highest priority on that position
        string[i] = str(position_pseudocode)
        # increment the position
        i += 1

    # join the list of characters into a string
    final_string = ''.join(string)
    return final_string




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


read_name = None  # variable which will store the name of the read currently being processed


for line in sys.stdin:  # reading lines of standard iput
    if line == '\n': # if it is an empty line
        continue # skip to the next one

    # splitting list into separate columns
    items = line.split()
    # changing name structure
    items[5] = re.sub('__.*', '', items[5])


    # if read_name is not defied then the first line is being processed
    if not read_name:
        r_length = int(items[1])
        # create a new dictionary
        d = new_dictionary(dictionary_upper, dictionary_lower, r_length)
        # increment the new dictionary
        d = increment(d, items[2], items[3], items[9], items[5], dictionary_upper, dictionary_lower)
        # reassign the read_name variable and contiue to the next line
        read_name = items[0]
        continue


    # if the first line has already been processed
    else:
        # if the same read is being processed
        if read_name == items[0]:
            # just increment the dictionary
            d = increment(d, items[2], items[3], items[9], items[5], dictionary_upper, dictionary_lower)
            # and reassign the read_name
            read_name = items[0]
            r_length = int(items[1])
            # continue to teh next read
            continue

        # if the loop moved on to teh next read
        else:
            # filter out the current dictionary
            d = filter(d, dictionary_lengths)
            # create the pseudocoded string
            s = pseudocode_string(d, dictionary_priorities, r_length)
            # write the pseudocoded string to a file
            print('>'+read_name+'\n'+s)
            #output_file.write('>'+read_name+'\n'+s+'\n')
            # reassign the read_name variable
            read_name = items[0]
            # and rad length variable
            r_length = int(items[1])
            # craete new dictionary for this read
            d = new_dictionary( dictionary_upper, dictionary_lower, r_length)
            # increment the dictionary
            d = increment(d, items[2], items[3], items[9], items[5], dictionary_upper, dictionary_lower)


# because the last read will come out of the loop without being processed it has to be processed additionally out
# of the loop
# filtering the dictionary
d = filter(d, dictionary_lengths)
# creating the pseudocode
s = pseudocode_string(d, dictionary_priorities, r_length)
# printing the pseudocoded reads as standard output
print('>'+read_name+'\n'+s)
#output_file.write('>'+read_name+'\n'+s+'\n')

# close output file
#output_file.close()



