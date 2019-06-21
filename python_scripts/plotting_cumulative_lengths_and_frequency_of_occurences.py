#!/usr/bin/env python3

from optparse import OptionParser

# This script will create a tabular output which has four columns and number of rows coresponding
# to the number of intervals in the binning of data.
# Each row of the output will have a sum of array lengths present within that bin.

# Format of input data looks like:

# satellite ID  # array length  # pseudocode	# read length	# I/T   # read ID
# LasTR1_6_13	803		        D		        82931		    I	    LAS180910g231_r34309_ch371_1D
# LasTR4        427366		    p		        135009		    T	    LAS180910g231_r7729_ch440_1D
# LasTR5	    415		        r		        42597		    I	    LAS190221_r714_ch107_1D

# defining the arguments which can be passed to the script
arguments = OptionParser()

arguments.add_option('-i', '--in', dest='table', help='length table file')
arguments.add_option('-n', '--number', dest='number_of_bins', help='number of bins')
arguments.add_option('-s', '--size', dest='size_of_bins', help='size of the bins')
arguments.add_option('-o', '--out', dest='output_name', help='name of the output file')

(options, args) = arguments.parse_args()
if (options.table is None or options.number_of_bins is None or
        options.size_of_bins is None or options.output_name is None):  # if one of the arguments is not provided
        print('\n----------> A mandatory option is missing !\n')   # raise an error

        arguments.print_help()  # and provide the help
        exit(-1)  # exit the script

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# if there will be any lengths above the limiting_bin value,
# they will all be pooled into the last bin
limiting_bin = int(options.size_of_bins) * int(options.number_of_bins)

# sequence of maximum values of each separate bin
bins = []

for i in range(int(options.size_of_bins), limiting_bin+int(options.size_of_bins), int(options.size_of_bins)):
        bins.append(i)

# these two cumulative dictionaries will hold the sum of array lengths within bins for intact and truncated
# satellite arrays respectively
# each satellite group has it's own key and a list with the summed values
cumulative_intact = {}
cumulative_truncated = {}
# these two frequency dictionaries will hold the frequency counts of arrays within bins for intact and truncated
# satellite arrays respectively
frequency_intact = {}
frequency_truncated = {}
# each of the satellite groups will have a separate tabular output

# opening the tabular input as f
with open(options.table) as f:
        # iterating over all lines in f
        for line in f:
                # splitting the columns within the line into separate items
                items = line.split()

                # if the satellite array is classified as intact
                if items[4] == 'I':

                        # if the satellite is not present within the dictionary keys
                        if items[0] not in cumulative_intact.keys():

                                # create the dictionary entry in all four dictionaries
                                # where the key will be the satellite name, and the value is a list
                                # filled with 0, the length of the number of bins
                                cumulative_intact[items[0]] = [0] * len(bins)
                                frequency_intact[items[0]] = [0] * len(bins)
                                cumulative_truncated[items[0]] = [0] * len(bins)
                                frequency_truncated[items[0]] = [0] * len(bins)

                                # to find the index of the bin to which the array length belongs to
                                # divide the length with the bin size
                                bin_index = int(items[1]) / int(options.size_of_bins)

                                # if the index is a round number, it means that the array length is the same as
                                # the maximum value of a bin, so in order to place the array within the bin
                                # 0.1 is subtracted from the array length
                                if bin_index == int(bin_index):
                                        # the bin_index is reassigned
                                         bin_index = (int(items[1])-0.1) / int(options.size_of_bins)

                                # if the bin_index is larger than the length of the bins
                                # it means the array length is larger than the limiting bin
                                if bin_index > int(options.number_of_bins):

                                        # therefore the array is placed within the last bin
                                         cumulative_intact[items[0]][-1] = cumulative_intact[items[0]][-1] + int(items[1])
                                         frequency_intact[items[0]][-1] = frequency_intact[items[0]][-1] + 1

                                # if the bin_index is smaller than the number of bins
                                if bin_index < int(options.number_of_bins):

                                        # use the bin_index as an index for the list belonging to the corresponding
                                        # satellite group within the two dictionaries
                                        cumulative_intact[items[0]][int(bin_index)] = cumulative_intact[items[0]][int(bin_index)] + int(items[1])
                                        frequency_intact[items[0]][int(bin_index)] = frequency_intact[items[0]][int(bin_index)] + 1

                        # if the satellite is present within the dictionary keys
                        else:
                                # the steps described previously in case where the satellite is not present within the
                                # keys of the dictionary are repeated
                                # but without creating a new entry

                                bin_index = int(items[1]) / int(options.size_of_bins)

                                if bin_index == int(bin_index):
                                        bin_index = (int(items[1]) - 0.1) / int(options.size_of_bins)

                                if bin_index > int(options.number_of_bins):
                                        cumulative_intact[items[0]][-1] = cumulative_intact[items[0]][-1] + int(items[1])
                                        frequency_intact[items[0]][-1] = frequency_intact[items[0]][-1] + 1

                                if bin_index < int(options.number_of_bins):
                                        cumulative_intact[items[0]][int(bin_index)] = cumulative_intact[items[0]][int(bin_index)] + int(items[1])
                                        frequency_intact[items[0]][int(bin_index)] = frequency_intact[items[0]][int(bin_index)] + 1

                # if the satellite arrays is classified as trunacted
                else:
                        # repeat all the steps described previously in cases where the satellite is and is not
                        # within the keys of the dictionary

                        if items[0] not in cumulative_truncated.keys():

                                cumulative_truncated[items[0]] = [0] * len(bins)
                                frequency_truncated[items[0]] = [0] * len(bins)
                                cumulative_intact[items[0]] = [0] * len(bins)
                                frequency_intact[items[0]] = [0] * len(bins)

                                bin_index = int(items[1]) / int(options.size_of_bins)

                                if bin_index == int(bin_index):
                                        bin_index = (int(items[1])) - 0.1 / int(options.size_of_bins)

                                if bin_index > int(options.number_of_bins):
                                        cumulative_truncated[items[0]][-1] = cumulative_truncated[items[0]][-1] + int(items[1])
                                        frequency_truncated[items[0]][-1] = frequency_truncated[items[0]][-1] + 1

                                if bin_index < int(options.number_of_bins):
                                        cumulative_truncated[items[0]][int(bin_index)] = cumulative_truncated[items[0]][int(bin_index)] + int(items[1])
                                        frequency_truncated[items[0]][int(bin_index)] = frequency_truncated[items[0]][int(bin_index)] + 1

                        else:
                                bin_index = int(items[1]) / int(options.size_of_bins)

                                if bin_index == int(bin_index):
                                        bin_index = (int(items[1]) - 0.1) / int(options.size_of_bins)

                                if bin_index > int(options.number_of_bins):
                                        cumulative_truncated[items[0]][-1] = cumulative_truncated[items[0]][-1] + int(
                                                items[1])
                                        frequency_truncated[items[0]][-1] = frequency_truncated[items[0]][-1] + 1

                                if bin_index < int(options.number_of_bins):
                                        cumulative_truncated[items[0]][int(bin_index)] = cumulative_truncated[items[0]][int(bin_index)] + int(items[1])
                                        frequency_truncated[items[0]][int(bin_index)] = frequency_truncated[items[0]][int(bin_index)] + 1

# check if the sets of keys from the dictionaries match in the truncated and intact dictionaries
if set(cumulative_intact) == set(cumulative_truncated) and set(frequency_intact) == set(frequency_truncated):

        # iterate over keys
        for key in cumulative_intact.keys():
                # because the keys are the same in all dictionaries the iteration can be used for all four dictionaries

                # creating a new output name for each satellite
                output_name = options.output_name+'_'+key
                # openning the output file
                out = open(output_name, 'w')

                # the first row of each output is the same
                first_output_row = 'Intact_cumulative   Truncated_cumulative    Frequency_intact        Frequency_truncated\n'
                out.write(first_output_row)

                # write the elements of the lists fromthe four dictionaries in rows corresponding to bins
                for (i,t,m,n) in zip(cumulative_intact[key], cumulative_truncated[key], frequency_intact[key], frequency_truncated[key]):
                        out.write(str(i)+'\t'+str(t)+'\t'+str(m)+'\t'+str(n)+'\n')
