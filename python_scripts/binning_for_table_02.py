#!/usr/bin/env python3

from optparse import OptionParser

# defining the arguments which can be passed to the script
arguments = OptionParser()

arguments.add_option('-n', '--number', dest='number', help='number of bins')
arguments.add_option('-s', '--size', dest='size', help='size of the bins')
arguments.add_option('-t', '--table', dest='table', help='table with numbers to bin')
arguments.add_option('-g', '--groups',  dest='groups', help='list of groups in format: -g a,b,c')
arguments.add_option('-o', '--out',  dest='output', help='name of output file')

options, args = arguments.parse_args()

if not options.number or not options.size or not options.table or not options.output or not options.groups:  # if one of the arguments is not provided
        print('\n----------> A mandatory option is missing !\n')   # raise an error

        arguments.print_help()  # and provide the help
        exit(-1)  # exit the script

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script accepts one additional argument: --groups.
# Each group is a separate satellite, and will have a separate output file.
# Format of input data looks like:

# satellite ID  # array length  # pseudocode	# read length	# I/T   # read ID
# LasTR1_6_13	803		D		82931		I	LAS180910g231_r34309_ch371_1D
# LasTR4	27366		p		135009		T	LAS180910g231_r7729_ch440_1D
# LasTR5	415		r		42597		I	LAS190221_r714_ch107_1D

# This script will eliminate duplicate reads if found within one group.

# the satellites variable is a list with the satellite IDs over which the script will iterate
satellites = options.groups.split(',')

# iterating over each satellite in the satellites list
for sat in satellites:

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

    read_name = []

    # opening the table file
    with open(options.table) as t:
    # iterating over each line in file
        for line in t:
	        # the columns in the table file are split into a list with separate items
            items = line.split()

	        # if the satellite ID from the file does not match the one being iterated over
            if items[0] != sat:
                # continue to next line
                continue

	        # if it does match
            else:

		        # check if the read ID already exists in the read_name list
                if items[5] in read_name:
		        # if it does, continue to next line
                    continue

		        # if it does not
                else:
		        # add the read ID to the read_name list
                    read_name.append(items[5])

                    # if the length of the read is larger than or equal to the limiting bin
                    if int(items[3]) >= limiting_bin:

                        # add the read length to the last bin
                        binning_data[limiting_bin].append(int(items[3]))

  		    # if it is not
                    else:

                        # i will be the variable used for counting, set it to 0
                        i = 0
                        # divide read length with a value from intervals[i]
			            # while the result of this division is smaller than 1
                        while int(items[3])/intervals[i] > 1:
                            # increment the i variable by 1
                            i += 1

                        # when the value of teh division is 1 or larger
                        # add the length of read to the corresponding bin
                        binning_data[intervals[i]].append(int(items[3]))

    # opening the output file
    out = open(options.output + '_' + sat, 'w')
    # iterating over keys in the sorted binning_data dictionary
    for key in sorted(binning_data.keys()):
        # into the output file, write: key and the sum of all lengths of that bin
        out.write(str(key) + '\t' + str(sum(binning_data[key])) + '\n')

