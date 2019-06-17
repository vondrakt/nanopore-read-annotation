#!/usr/bin/env python3

from optparse import OptionParser

# defining the arguments which can be passed to the script
arguments = OptionParser()

arguments.add_option('-n', '--number', dest='number', help='number of bins')
arguments.add_option('-s', '--size', dest='size', help='size of the bins')
arguments.add_option('-t', '--table', dest='table', help='table with numbers to bin')
arguments.add_option('-c', '--column', dest='column', help='column to use for binning, COUNTING FROM 0!')
arguments.add_option('-o', '--out',  dest='output', help='name of output file')

options, args = arguments.parse_args()

if not options.number or not options.size or not options.table or not options.output or not options.column:  # if one of the arguments is not provided
        print('\n----------> A mandatory option is missing !\n')   # raise an error

        arguments.print_help()  # and provide the help
        exit(-1)  # exit the script

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script works for files with multiple columns.
# Using the argument --column it is possible to choose which column will be binned.

# I have used it for data with a structure:

# Read ID                       Read length
# LAS180910g231_r2333_ch115_1D	42096
# LAS180910g231_r2462_ch439_1D	30780
# LAS180910g231_r2492_ch115_1D	92780

# If there are duplicate reads within the files,
# they will not be filtered out.
print('This script does not take into consideration if there are duplicate sequences.')

# if there will be any lengths above the limiting_bin value,
# they will all be pooled into the last bin
limiting_bin = int(options.number) * int(options.size)

# data frame that contains the maximum values of each separate bin as keys
# and lists filled with 0 as values
binning_data = {}
# sequence of maximum values of each separate bin
intervals = []

for i in range(int(options.size), limiting_bin + int(options.size), int(options.size)):
    binning_data[i] = [0]
    intervals.append(i)

# the index variable stores which column to be used for binning
index = int(options.column)

# opening the table file
with open(options.table) as t:
    # iterating over each line of file
    for line in t:

        # the columns in the table file are split into a list with separate items
        items = line.split()
        # if the read length is larger than or equal to the limiting_bing
        if int(items[index]) >= limiting_bin:
            # add the value to the largest bin in the dataframe
            binning_data[limiting_bin].append(int(items[index]))

        # if it is not
        else:
            # i will be the variable used for counting, set it to 0
            i = 0
            # divide read length with a value from intervals[i]
            # while the result of this division is smaller than 1
            while int(items[index])/intervals[i] > 1:
                i += 1

            # when the value of teh division is 1 or larger
            # add the length of read to the corresponding bin
            binning_data[intervals[i]].append(int(items[index]))

# opening the output file
out = open(options.output, 'w')

# iterating over keys in the sorted binning_data dictionary
for key in sorted(binning_data.keys()):
    # into the output file, write: key and the sum of all lengths of that bin
    out.write(str(key) + '\t' + str(sum(binning_data[key])) + '\n')
