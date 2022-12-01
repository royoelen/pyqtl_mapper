
import numpy

class Covariates:


    def __init__(self, covariate_file_location, covariates_to_use=None):
        #self.covariate_to_index = None
        self.donor_to_index = None
        self.covariates = None
        self.set_up_covariates(covariate_file_location, covariates_to_use)

    def set_up_covariates(self, covariate_file_location, covariates_to_use=None, donor_column='donor'):
        # we will initialize the mappings for where the donors and covariates will be in the 2d numpy array
        self.covariates = {}
        self.donor_to_index = {}
        # we need to know if we are at the first line
        at_header = True
        # we need to keep an index of where each donor is in the numpy array
        index_donor = 0
        # we using the header to keep track of where in our line is which covariate
        column_mapping = {}
        # and we have the donor column index
        donor_column_index = None
        # open the covariate file
        with open(covariate_file_location) as file_connection:
            for line in file_connection:
                # we'll split the header
                line_split = line.split()
                # at the header, we need to set up some stuff
                if at_header:
                    # will set up the mappings for which column is which covariate
                    for column_index in range(0, len(line_split), 1):
                        # get the column name
                        variable_name = line_split[column_index]
                        # depending on if we want that covariate, we will put it in our column mapping, if no covariates to use are supplied, we just do all
                        if (covariates_to_use is None or variable_name in covariates_to_use) and variable_name != donor_column:
                            # get the index of the column where the variable is
                            column_mapping[variable_name] = column_index
                            # create a numpy array for the covariate
                            array_covar = numpy.empty(1)
                            # add that to the dictionary that has an array for each covariate, keyed by the covariate name
                            self.covariates[variable_name] = array_covar
                        elif variable_name == donor_column:
                            # get the index of the column where the variable is
                            donor_column_index = column_index

                    # okay, that is enough work for the header
                    at_header = False
                # after the header we start adding data
                else:
                    # get the donor
                    donor_name = line_split[donor_column_index]
                    # map the donor to the index
                    self.donor_to_index[donor_name] = index_donor
                    # check each covariate that we want to use
                    for covariate, index in column_mapping.items():
                        # check if the array is still large enough
                        if index_donor >= self.covariates[covariate].size:
                            # increase the size if neccessary
                            self.covariates[covariate] = numpy.resize(self.covariates[covariate], 10 * self.covariates[covariate].size)
                        # add to that array
                        self.covariates[covariate][index_donor] = line_split[index]
                    # a new donor is on each line, so increase the index
                    index_donor = index_donor + 1

    def get_valid_donors(self):
        # TODO actually implement
        print('get_valid_donors not implemented')
        return list(self.donor_to_index.keys())