
import QTLConfig
import ProbeLocation
import SNPLocation

from numpy.core.fromnumeric import ndim
from sklearn import linear_model
from scipy import stats
import numpy as numpy
import pandas as pandas
from multiprocessing import Process, Queue


class QTLMapper:

    def __init__(self,  snp_file_location, probe_file_location, covariates_file_location=None, covariates_to_use=None, snp_positions_file_location = None, probe_positions_file_location=None, use_model='linear', confinements_snp_location=None, confinements_probe_location=None, confinements_snp_probe_pairs_location=None, maf=0.01):
        self.qtl_config = QTLConfig.QTLConfig()
        self.qtl_config.snp_file_location = snp_file_location
        self.qtl_config.probe_file_location = probe_file_location
        self.qtl_config.covariates_file_location = covariates_file_location
        self.qtl_config.covariates_to_use = covariates_to_use
        self.qtl_config.snp_positions_file_location = snp_positions_file_location
        self.qtl_config.probe_positions_file_location = probe_positions_file_location
        self.qtl_config.use_model = use_model
        self.qtl_config.confinements_snp_location = confinements_snp_location
        self.qtl_config.confinements_probe_location = confinements_probe_location
        self.qtl_config.confinements_snp_probe_pairs_location = confinements_snp_probe_pairs_location
        self.qtl_config.maf = maf

        # initialize variables of class
        self.snp_confinement = None
        self.probe_confinement = None
        self.snp_probe_confinement = None
        self.probe_locations = None
        self.snp_locations = None
        self.covariates = None

        # connection to file
        self.active_file_connection = None

        # set confinement files
        self.set_confinements()

        # set the location files
        self.set_locations()

        # set the covariates
        self.set_covariates()


    def set_confinements(self):
        # read the confinement file for SNPs
        if self.qtl_config.confinements_snp_location is not None:
            self.snp_confinement = pandas.read_csv(self.qtl_config.confinements_snp_location, sep = '\t')
        # read the confinement file for probes
        if self.qtl_config.confinements_probe_location is not None:
            self.probe_confinement = pandas.read_csv(self.qtl_config.confinements_probe_location, sep = '\t')
        # read the confinement file for snp-probe combinations
        if self.qtl_config.confinements_snp_probe_pairs_location is not None:
            self.snp_probe_confinement = pandas.read_csv(self.qtl_config.confinements_snp_probe_pairs_location, sep = '\t')


    def set_locations(self):
        # set the locations of the position files
        if self.qtl_config.snp_positions_file_location is not None:
            self.snp_locations = SNPLocation.SNPLocation(self.qtl_config.snp_positions_file_location)
        if self.qtl_config.probe_positions_file_location is not None:
            self.probe_locations = ProbeLocation.ProbeLocation(self.qtl_config.probe_positions_file_location)


    def set_covariates(self):
        print('would set covariates')

    def check_maf(self, genotypes):
        # doing 0,1,2
        allele_frequency_alt = sum(genotypes) / (len(genotypes) * 2)
        minor_allele_frequency = allele_frequency_alt
        if minor_allele_frequency > 0.5:
            minor_allele_frequency = 1 - minor_allele_frequency
        # check against the set MAF
        if minor_allele_frequency > self.qtl_config.maf:
            return True
        else:
            return False


    def perform_regression(self, X, y):
        # create a linear regression model
        model = LinearRegression()

        # do a fit
        model = model.fit(X, y)

        # coefficient of determination
        r_sq = model.score(X, y)
        # t value
        t = model.t
        # p value
        p = model.p


    def decide_mapping_snp(self, snp_id):
        # if there is no confinement, I guess we'll try
        if self.snp_confinement is None and self.snp_probe_confinement is None:
            return True
        # if we have just a SNP confinement, check if the SNP is in there
        elif self.snp_confinement is not None and self.snp_probe_confinement is None:
            # check if it is in the snp confinement
            if len(self.snp_confinement[self.snp_confinement[0].str.contains(snp_id)]) > 0:
                return True
            else:
                return False
        # if there is a snp-probe confinement, also check for the SNP
        elif self.snp_confinement is  None and self.snp_probe_confinement is not None:
            # check if it is in the snp confinement
            if len(self.snp_probe_confinement[self.snp_probe_confinement[0].str.contains(snp_id)]) > 0:
                return True
            else:
                return False
        # for both, I guess we will check both
        elif self.snp_confinement is not None and self.snp_probe_confinement is not None:
            # check if it is in the snp confinement and the snp-probe confinement
            if (len(self.snp_confinement[self.snp_confinement[0].str.contains(snp_id)]) > 0) and (len(self.snp_probe_confinement[self.snp_probe_confinement[0].str.contains(snp_id)]) > 0):
                return True
            else:
                return False

    def check_cis_distance(self, snp_id, probe_id):
        # if we don't have a cis distance, we'll just do anything
        if self.qtl_config.cis_dist is None:
            return True
        # zero means also means we'll just do anything
        elif self.qtl_config.cis_dist == 0:
            return True
        # apparently we need to check the cis distance
        else:
            # check if we have a positions file set
            if self.probe_locations is not None and self.snp_locations is not None:
                # get the snp data
                snp_data = self.snp_locations.get_snp_position(snp_id)
                # get the probe data
                probe_data = self.probe_locations.get_probe_position(probe_id)
            else:
                print('using cis distance, but not both snp and probe locations are supplied')


    def confinement_inclusion(self, snp_id, probe_id):
        # if there is no confinement, I guess we'll try
        if self.probe_confinement is None and self.snp_probe_confinement is None:
            return True
        # if there is a probe confinement, check if it is in there
        elif self.probe_confinement is not None and self.snp_probe_confinement is None:
            # now we have to see if the probe is in the inclusion list
            if(len(self.probe_confinement[self.probe_confinement[0].str.contains(snp_id)]) > 0):
                return True
            else:
                return False
        # if there is a SNP probe confinement, check if there is a match
        elif self.snp_probe_confinement is not None:
            return False


    def map_qtl(self, genotype_data, probe_data, genotype_metadata, probe_metadata):

        self.perform_regression(X, y)


    def do_mapping_snp(self, snp_id, donors_snps, genotypes, donor_offset=1):
        # open connection with the gene file, this we will do for every SNP, to be able to do it in parallel
        probe_active_file_connection = open(self.qtl_config.confinements_snp_location)
        # we will get the header
        donors_probes = None
        # and remember it was the header
        is_header = True
        # check each line
        for line in probe_active_file_connection:
            # split into pieces
            data_row = line.split(sep='\t')
            # turn into numpy
            data_row = numpy.array(data_row)
            # check if it is the header
            if is_header is not True:
                # the first entry is the identifier
                probe_id = data_row[0]
                # check for gene confinement
                if self.confinement_inclusion:
                    # next check for the cis distance
                    if self.check_cis_distance(snp_id, probe_id): # okay, now it is worth the effort to check the probe contents
                        # replace empty values with nans
                        data_row = numpy.char.replace(data_row, 'na', 'nan')
                        data_row = numpy.char.replace(data_row, '', 'nan')
                        # turn into floats
                        probes = data_row[1:].astype(numpy.float)
                        # get the valid probes
                        indices_valid_probes = numpy.where((numpy.isfinite(probes) == False) & (numpy.isnan(probes) == False))
                        # get the donors with these valid probes
                        valid_probe_donors = donors_probes[indices_valid_probes]
                        # overlap with the valid SNP donors
                        common_donors = numpy.intersect1d(donors_snps, valid_probe_donors, assume_unique=True)
                        # prepare sorted genotype and probe arrays
                        genotypes_sorted = numpy.empty(shape=len(common_donors))
                        probes_sorted = numpy.empty(shape=len(common_donors))
                        # by using indices
                        sorted_index = 0
                        # fill these
                        for donor in common_donors:
                            # find the index of the genotype
                            index_genotype = numpy.where(donors_snps == donor)
                            # find the index of the probe
                            index_probe = numpy.where(donors_probes == donor)
                            # get the genotype
                            geno_donor = genotypes[index_genotype]
                            # get the probe
                            probe_donor = probes[index_probe]
                            # add to the numpy array
                            genotypes_sorted[sorted_index] = geno_donor
                            probes_sorted[sorted_index] = probe_donor
                            # increase the index
                            sorted_index = sorted_index + 1
                        # now check if we are okay with the MAF
                        if self.check_maf(genotypes_sorted):
                            print('blie')
                        else:
                            print('blah')

            else:
                # need to set the header only once
                donors_probes = data_row[range(donor_offset, len(data_row, 1))]
                # subset the
                is_header = False


    def perform_mapping(self, donor_offset=1):
        # open connection with SNP file
        self.active_file_connection = open(self.qtl_config.snp_file_location)
        # we will get the header
        donors_snps = None

        # and remember it was the header
        is_header = True
        # check each line
        for line in self.active_file_connection:
            # split into pieces
            data_row = line.split(sep='\t')
            # turn into numpy
            data_row = numpy.array(data_row)
            # check if it is the header
            if is_header is not True:
                # the first entry is the identifier
                snp_id = data_row[0]
                # check if we want to map the SNP
                if self.decide_mapping_snp(snp_id):
                    # get the genotypes
                    genotypes_list = data_row[list(range(donor_offset, len(data_row), 1))] # there can be other info fields in addition to the snp id, with the offset you can take this in consideration
                    # to numpy for speed
                    genotypes = numpy.array(genotypes_list)
                    # replace the empty entries with -1
                    genotypes = numpy.char.replace(genotypes, '.', '-1')
                    genotypes = numpy.char.replace(genotypes, 'nan', '-1')
                    genotypes[genotypes == '']='-1'
                    # as floats, not strings of course
                    genotypes = genotypes.astype(numpy.float)
                    # check which donors have correct values in the SNPs
                    indices_valid_genotype = (genotypes >= 0).nonzero()
                    # subset the donor names and their genotypes
                    donors_valid_snp = numpy.take(donors_snps, indices_valid_genotype)
                    genotypes_valid_snp = numpy.take(genotypes, indices_valid_genotype)
                    # perform mapping with this SNP
                    self.do_mapping_snp(snp_id, donors_valid_snp, genotypes_valid_snp)



            else:
                # need to set the header only once
                donors_snps = data_row[list(range(donor_offset, len(data_row), 1))]
                is_header = False




    def read_chunks_once(self, file_connection, chunk_size=1024):
        '''

        :param file_connection:
        :param chunk_size:
        :return:
        '''
        while True:
            data = file_connection.read(chunk_size)
            if not data:
                break
            yield data



class LinearRegression(linear_model.LinearRegression):
    """
    LinearRegression class after sklearn's, but calculate t-statistics
    and p-values for model coefficients (betas).
    Additional attributes available after .fit()
    are `t` and `p` which are of the shape (y.shape[1], X.shape[1])
    which is (n_features, n_coefs)
    This class sets the intercept to 0 by default, since usually we include it
    in X.
    """

    def __init__(self, *args, **kwargs):
        if not "fit_intercept" in kwargs:
            kwargs['fit_intercept'] = False
        super(LinearRegression, self) \
            .__init__(*args, **kwargs)

    def fit(self, X, y, n_jobs=1):
        self = super(LinearRegression, self).fit(X, y, n_jobs)

        sse = numpy.sum((self.predict(X) - y) ** 2, axis=0) / float(X.shape[0] - X.shape[1])
        # standard error
        se = numpy.array([numpy.sqrt(numpy.diagonal(sse * numpy.linalg.inv(numpy.dot(X.T, X))))])

        self.t = self.coef_ / se
        self.p = 2 * (1 - stats.t.cdf(numpy.abs(self.t), y.shape[0] - X.shape[1]))
        return self
