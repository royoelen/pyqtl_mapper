import QTLConfig

import numpy as numpy

class ProbeLocation:


    def __init__(self, probe_locations_file_location):
        # we have a scheme where the genes are keys for location arrays
        self.probe_to_index = None
        # the there will be three numpy arrays containing the variables we need
        self.chromosome_loc = None
        self.probe_start_loc = None
        self.probe_stop_loc = None


    def set_up_locations(self, probe_locations_file_location):
        # create the index file
        self.probe_to_index = {}
        # we will use numpy arrays for the variables
        self.chromosome_loc = numpy.empty(1, dtype='S3')
        self.probe_start_loc = numpy.empty(1, dtype='uint64')
        self.probe_stop_loc = numpy.empty(1, dtype='uint64')

        # there is header
        at_header = True
        # and we will keep an index
        index = 0
        # this file may be large, so we will go line by line
        with open(probe_locations_file_location) as file_connection:
            for line in file_connection:
                # skip header
                if at_header:
                    at_header = False
                else:
                    # we grow the numpy arrays, depending on the number of items
                    if index >= self.chromosome_loc.size:
                        # increases sizes by ten
                        self.chromosome_loc = numpy.resize(self.chromosome_loc, 10 * self.chromosome_loc.size)
                        self.probe_start_loc = numpy.resize(self.probe_start_loc, 10 * self.probe_start_loc.size)
                        self.probe_stop_loc = numpy.resize(self.probe_stop_loc, 10 * self.probe_stop_loc.size)
                    # split the line into three
                    split_line = line.split()
                    # set the index for this probe
                    self.probe_to_index[split_line[0]] = index
                    self.chromosome_loc[index] = split_line[1]
                    self.probe_start_loc[index] = int(split_line[2])
                    self.probe_to_index[index] = int(split_line[3])

    def get_probe_position(self, probe_name):
        # get the index for the probe
        probe_index = self.probe_to_index[probe_name]
        # put result in dictionary
        probe_info = {
            'geneid' : probe_name,
            'chr' : self.chromosome_loc[probe_index],
            'left' : self.probe_start_loc[probe_index],
            'right' : self.probe_stop_loc[probe_index]
        }
        return probe_info