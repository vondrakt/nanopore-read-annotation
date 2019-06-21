#!/usr/bin/env python3

from optparse import OptionParser
import sys

# defining the arguments that need to be passed to the script
arguments = OptionParser()

arguments.add_option('-b', '--bit', dest='bit_score', help='minimum bit score to be filtered for')
arguments.add_option('-x', '--extension', dest='extension_percentage', help='maximum hit length based on the percentage of reference database sequence length')

(options, args) = arguments.parse_args()

if (options.bit_score is None or options.extension_percentage is None):
    # if one of the arguments is missing
    print('\n----------> A mandatory option is missing !\n')  # raise an error
    arguments.print_help()  # and provide the help
    exit(-1)  # exit the script

# read lines of the standard input
for line in sys.stdin:
    # separate the lines into individual items based on the tab	
    items = line.split()

    # store the 1.23x length of the reference sequence in the extension_length variable
    extension_length = int(items[6]) * float(options.extension_percentage)

    # if the length of the hit is smaller than the extension_length and bitcore of the hit larger than the given limiting value
    if int(items[3]) < extension_length and int(items[12]) > int(options.bit_score):
        # print the line
        print(line)
