import QTLConfig

import numpy as numpy

class SNPLocation:


    def __init__(self, snp_locations_file_location):
        # we have a scheme where the genes are keys for location arrays
        self.snp_to_index = None
        # the there will be two numpy arrays containing the variables we need
        self.chromosome_loc = None
        self.snp_loc = None



    def set_up_locations(self, snp_locations_file_location):
        # create the index file
        self.snp_to_index = {}
        # we will use numpy arrays for the variables
        self.chromosome_loc = numpy.empty(1, dtype='S3')
        self.snp_loc = numpy.empty(1, dtype='uint64')

        # there is header
        at_header = True
        # and we will keep an index
        index = 0
        # this file may be large, so we will go line by line
        with open(snp_locations_file_location) as file_connection:
            for line in file_connection:
                # skip header
                if at_header:
                    at_header = False
                else:
                    # we grow the numpy arrays, depending on the number of items
                    if index >= self.chromosome_loc.size:
                        # increases sizes by ten
                        self.chromosome_loc = numpy.resize(self.chromosome_loc, 10 * self.chromosome_loc.size)
                        self.snp_loc = numpy.resize(self.snp_loc, 10 * self.snp_loc.size)
                    # split the line into three
                    split_line = line.split()
                    # set the index for this probe
                    self.snp_to_index[split_line[0]] = index
                    self.chromosome_loc[index] = split_line[1]
                    self.snp_loc[index] = int(split_line[2])

    def get_snp_position(self, snp_name):
        # get the index for the probe
        snp_index = self.snp_to_index[snp_name]
        # put result in dictionary
        snp_info = {
            'snpid' : snp_name,
            'chr' : self.chromosome_loc[snp_index],
            'pos' : self.snp_loc[snp_index],
        }
        return snp_info
