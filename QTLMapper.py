
import QTLConfig
import ProbeLocation
import SNPLocation
import Covariates

from numpy.core.fromnumeric import ndim
from sklearn import linear_model
from scipy import stats
import numpy as numpy
import pandas as pandas
import statsmodels.stats.multitest as multitest
import warnings
from multiprocessing import Process, Queue, Pool


class QTLMapper:
    """
    Class to perform QTL mapping
    """

    def __init__(self,  snp_file_location, probe_file_location, output_location, covariates_file_location=None, covariates_to_use=None, snp_positions_file_location = None, probe_positions_file_location=None, use_model='linear', confinements_snp_location=None, confinements_probe_location=None, confinements_snp_probe_pairs_location=None, maf=0.01, cis_dist=1000000):
        """

        :param snp_file_location: location of the SNP file, expected to be tsv with variant name, then variants per donor
        :param probe_file_location: location of the probe file, expected to tsv with probe name, then value per donor
        :param output_location: output tsv location of the resulting mapping
        :param covariates_file_location: location of the covariates file, expected to be tsv
        :param covariates_to_use: covariates to use, need to match the column names in covariates file
        :param snp_positions_file_location: position annotation file for variants
        :param probe_positions_file_location: position annotation file for probes
        :param use_model: the type of regression model to use
        :param confinements_snp_location: optional SNP confinement, variant names need to match the snp file variant names
        :param confinements_probe_location: optional probe confinement, probe names need to match the probe file probe names
        :param confinements_snp_probe_pairs_location: optional SNP-probe confinement, , variant names need to match the snp file variant names, probe names need to match the probe file probe names
        :param maf: minimal minor allele frequency cutoff to consider a variant worth performing QTL mappinf for
        :param cis_dist: max distance from the probe start or stop, that a SNP is still considered within cis distance for
        """
        self.qtl_config = QTLConfig.QTLConfig()
        self.qtl_config.snp_file_location = snp_file_location
        self.qtl_config.probe_file_location = probe_file_location
        self.qtl_config.output_location = output_location
        self.qtl_config.covariates_file_location = covariates_file_location
        self.qtl_config.covariates_to_use = covariates_to_use
        self.qtl_config.snp_positions_file_location = snp_positions_file_location
        self.qtl_config.probe_positions_file_location = probe_positions_file_location
        self.qtl_config.use_model = use_model
        self.qtl_config.confinements_snp_location = confinements_snp_location
        self.qtl_config.confinements_probe_location = confinements_probe_location
        self.qtl_config.confinements_snp_probe_pairs_location = confinements_snp_probe_pairs_location
        self.qtl_config.maf = maf
        self.qtl_config.cis_dist = cis_dist

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
        """
        set the confinements from the input
        :return: None
        """
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
        """
        set the locations from the input
        :return: None
        """
        # set the locations of the position files
        if self.qtl_config.snp_positions_file_location is not None:
            self.snp_locations = SNPLocation.SNPLocation(self.qtl_config.snp_positions_file_location)
        if self.qtl_config.probe_positions_file_location is not None:
            self.probe_locations = ProbeLocation.ProbeLocation(self.qtl_config.probe_positions_file_location)

    def set_covariates(self):
        """
        set the covariates from the input
        :return: None
        """
        if self.qtl_config.covariates_file_location is not None:
            self.covariates = Covariates.Covariates(self.qtl_config.covariates_file_location, self.qtl_config.covariates_to_use)

    def check_maf(self, genotypes):
        """
        check the minor allele frequency of given genotypes
        :param genotypes: genotypes in 0 homozygous wildtype, 1 heterozygote, 2 homozygous mutant format
        :return: whether the minor allele frequency was below the cutoff
        """
        # doing 0,1,2
        allele_frequency_alt = sum(genotypes) / (len(genotypes) * 2)
        minor_allele_frequency = allele_frequency_alt
        # minor allele 0.1 is the same as 0.9 actually, so we need to check both ways
        if minor_allele_frequency > 0.5:
            minor_allele_frequency = 1 - minor_allele_frequency
        # check against the set MAF
        if minor_allele_frequency > self.qtl_config.maf:
            return True
        else:
            return False

    def perform_regression(self, X, y):
        """
        perform regression
        :param X: a pandas dataframe containing the variates to use for the prediction
        :param y: the predicted value
        :return: a LinearRegression object, resulting from running the regression
        """
        # create a linear regression model
        model = LinearRegression()
        # do a fit
        model = model.fit(X, y)
        return model

    def decide_mapping_snp(self, snp_id):
        """
        decide on whether or not consider the supplied SNP for QTL mapping
        :param snp_id: the variant ID of the supplied SNP
        :return: whether or not the supplied SNP should be considered for QTL mapping
        """
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
        """
        check whether or not a given variant is within cis distance of a given probe
        :param snp_id: the variant ID of the SNP
        :param probe_id: the name of the probe
        :return: whether or not the given SNP was within cis distance of the given probe
        """
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
                # check if there was data for the probe
                if probe_data is not None and snp_data is not None:
                    # check if within cis distance
                    max_loc_left = None
                    max_loc_right = None
                    # depending on the strand, left and right can be different
                    if probe_data['right'] > probe_data['left']:
                        max_loc_left = probe_data['left'] - self.qtl_config.cis_dist
                        max_loc_right = probe_data['right'] + self.qtl_config.cis_dist
                    elif probe_data['right'] < probe_data['left']:
                        max_loc_left = probe_data['left'] + self.qtl_config.cis_dist
                        max_loc_right = probe_data['right'] - self.qtl_config.cis_dist
                    # check if the probe is somewhere between these locations, on the same chromosome
                    #print(' '.join([str(probe_data['chr']), str(snp_data['chr']), str(snp_data['pos']), str(max_loc_left), str(max_loc_right)]))
                    if probe_data['chr'] == snp_data['chr'] and max_loc_left <= snp_data['pos'] <= max_loc_right:
                        return True
                    else:
                        return False

                elif probe_data is None and snp_data is not None:
                    warnings.warn(' '.join(['probe location not found in location file, skipping entry', str(snp_id), ';', str(snp_data), ',', str(probe_id), ':', str(probe_data)]))
                    return False
                elif probe_data is not None and snp_data is None:
                    warnings.warn(' '.join(['snp location not found in location file, skipping entry', str(snp_id), ';', str(snp_data), ',', str(probe_id), ':', str(probe_data)]))
                    return False
                else:
                    # should never happen
                    #warnings.warn(' '.join(['probe and/or SNP location not found in location file, skipping entry', str(snp_id), ';', str(snp_data), ',', str(probe_id), ':', str(probe_data)]))
                    return False
            else:
                print('using cis distance, but not both snp and probe locations are supplied, this will cause all probes to be skipped')
                return False

    def confinement_inclusion(self, snp_id, probe_id):
        """
        check if a snp-probe combination is present in our confinement. If the confinement files were not supplied, they are not used, and TRUE will be returned

        :param snp_id: the variant identifier
        :param probe_id: the probe identifier
        :return: whether the combination was present in the confinement. If the confinement files were not supplied, they are not used, and TRUE will be returned
        """
        # if there is no confinement, I guess we'll try
        if self.probe_confinement is None and self.snp_probe_confinement is None and self.snp_confinement is None:
            return True
        # if there is only a probe confinement
        elif self.probe_confinement is not None and self.snp_probe_confinement is None and self.snp_confinement is None:
            # now we have to see if the probe is in the inclusion list
            if(len(self.probe_confinement[self.probe_confinement[0].str.contains(probe_id)]) > 0):
                return True
            else:
                return False
        # if there is only a SNP confinement
        elif self.probe_confinement is None and self.snp_probe_confinement is None and self.snp_confinement is not None:
            # now we have to see if the probe is in the inclusion list
            if(len(self.snp_confinement[self.snp_confinement[0].str.contains(snp_id)]) > 0):
                return True
            else:
                return False
        # if there is a only a snp-probe confinement, check if it is in there
        elif self.probe_confinement is None and self.snp_probe_confinement is not None and self.snp_confinement is None:
            # now we have to see if the snp-probe combination is in the inclusion list
            if self.snp_probe_confinement[(self.snp_probe_confinement['snp'] == snp_id & self.snp_probe_confinement['probe'] == probe_id)].shape[1] > 0:
                return True
            else:
                return False
        # if there is a snp-probe confinement, and a snp confinement
        elif self.probe_confinement is None and self.snp_probe_confinement is not None and self.snp_confinement is not None:
            if self.snp_probe_confinement[(self.snp_probe_confinement['snp'] == snp_id & self.snp_probe_confinement['probe'] == probe_id)].shape[1] > 0 and len(self.snp_confinement[self.snp_confinement[0].str.contains(snp_id)]) > 0:
                return True
            else:
                return False
        # if there is a snp-probe confinement, and a probe confinement
        elif self.probe_confinement is not None and self.snp_probe_confinement is not None and self.snp_confinement is None:
            if self.snp_probe_confinement[(self.snp_probe_confinement['snp'] == snp_id & self.snp_probe_confinement['probe'] == probe_id)].shape[1] > 0 and len(self.probe_confinement[self.probe_confinement[0].str.contains(probe_id)]) > 0:
                return True
            else:
                return False
        # if all three confinements are present
        elif self.probe_confinement is not None and self.snp_probe_confinement is not None and self.snp_confinement is not None:
            if self.snp_probe_confinement[(self.snp_probe_confinement['snp'] == snp_id & self.snp_probe_confinement['probe'] == probe_id)].shape[1] > 0 and len(self.snp_confinement[self.snp_confinement[0].str.contains(snp_id)]) > 0 and len(self.probe_confinement[self.probe_confinement[0].str.contains(probe_id)]) > 0:
                return True
            else:
                return False

    def map_qtl(self, probe_values, covariate_table, genotypes_sorted):
        """
        perform the QTL mapping

        :param probe_values: the values for the probe, sorted as the covariate table and genotypes
        :param covariate_table: the covariate table, sorted as the probe values and genotypes
        :param genotypes_sorted: the genotypes, sorted as the covariate table and probes
        :return: a dictionary containing the result of the mapping
        """
        # initilialize models
        model_null = None
        model_genotypes = None
        # we can only do this if there is metadata
        if covariate_table is not None:
            # perform regression without the genotypes
            model_null = self.perform_regression(covariate_table, probe_values)
            # perform the regression with the genotypes
            covariate_table['snp'] = genotypes_sorted
        else:
            covariate_table = pandas.DataFrame({'snp':genotypes_sorted})
        # regresss with genotypes
        model_genotypes = self.perform_regression(covariate_table, probe_values)
        # extract information from the genotype model

        # coefficient of determination
        r_sq = model_genotypes.score(covariate_table, probe_values)
        # t value
        t = model_genotypes.t[0][0]
        # p value
        p = model_genotypes.p[0][0]

        result = {'r_sq' : r_sq, 't' : t, 'p' : p}

        # now if we did a null regression, we can see how well we predict with the genotypes and without the genotypes
        if model_null is not None:
            rss_null = model_null.ssr
            rss_geno = model_genotypes.ssr
            # TODO further implementation

        return result

    def get_genotype_metadata_as_pandas(self, donors):
        """
        get the metadata for the given donors, in that donor order
        :param donors: the donors to fetch the metadata for
        :return: the metadata for the given donors in that donor order
        """
        # initialize value
        covars_pandas = None
        # we can only fill in the value if it is not none
        if self.covariates is not None:
            # turn the metadata into a pandas dataframe
            covars_pandas = pandas.DataFrame({covariate:list(np_array) for covariate,np_array in self.covariates.covariates.items()})
            # get the indices of the donors using list comprehension
            covariate_donor_indices = [self.covariates.donor_to_index.get(donor) for donor in donors]
            # subset to the indices of the donors
            covars_pandas = covars_pandas.iloc(covariate_donor_indices)

        return covars_pandas


    def get_probe_metadata_as_pandas(self, donors):
        print('unimplemented')

    def do_mapping_snp(self, snp_id, donors_snps, genotypes, donor_offset=1):
        """
        perform the QTL mapping for a given variant

        :param snp_id: the variant ID
        :param donors_snps: the donors that we have the variant for
        :param genotypes: the genotypes of the donors
        :param donor_offset: the offset fo the genotype file where the genotypes start
        :return: list of dictionaries containing the result of each mapping
        """
        # open connection with the gene file, this we will do for every SNP, to be able to do it in parallel
        probe_active_file_connection = open(self.qtl_config.probe_file_location)
        # we will get the header
        donors_probes = None
        # and remember it was the header
        is_header = True
        # save the results
        results = []
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
                        # turn into floats
                        probes = data_row[1:].astype(numpy.float)
                        # get the valid probes
                        indices_valid_probes = numpy.where(numpy.isfinite(probes))
                        # get the donors with these valid probes
                        valid_probe_donors = donors_probes[indices_valid_probes]
                        # overlap with the valid SNP donors
                        common_donors = numpy.intersect1d(donors_snps, valid_probe_donors, assume_unique=True)
                        # also check with the covariates if we use those
                        if self.covariates is not None:
                            # get valid covariate donors
                            valid_covariate_donors = self.covariates.get_valid_donors()
                            # convert to numpy array
                            valid_covariate_donors = numpy.array(valid_covariate_donors)
                            # then go to common donors
                            common_donors = numpy.intersect1d(common_donors, valid_covariate_donors)
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
                            # get a Pandas dataframe with the metadata
                            total_covariates = self.get_genotype_metadata_as_pandas(common_donors)
                            # perform the regression
                            result_snp_gene = self.map_qtl(probes_sorted, total_covariates, genotypes_sorted)
                            # add the snp and gene as keys
                            result_snp_gene['snp'] = snp_id
                            result_snp_gene['probe'] = probe_id
                            # add to the results
                            results.append(result_snp_gene)
                        else:
                            warnings.warn(' '.join(['skipping due to MAF', snp_id]))
                    # if we are not in the cis distance
                    else:
                        pass
                # if the gene was no in the confinement
                else:
                    pass
            # this is the header
            else:
                # need to set the header only once
                donors_probes = data_row[range(donor_offset, len(data_row), 1)]
                # subset the
                is_header = False
        # return list of results
        return results

    def filter_and_map_snp(self, data_row, donor_offset, donors_snps, snp_id):
        """
        perform QTL mapping for a given SNP

        :param data_row: numpy array containing the genotypes of the donors
        :param donor_offset: there can be other info fields in addition to the snp id, with the offset you can take this in consideration
        :param donors_snps: the donors the genotypes are for
        :param snp_id: the variant ID under examination
        :return: list of dictionaries containing the results
        """
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
        result = self.do_mapping_snp(snp_id, donors_valid_snp, genotypes_valid_snp)
        return result

    def perform_mapping(self, donor_offset=1):
        # open connection with SNP file
        self.active_file_connection = open(self.qtl_config.snp_file_location)
        # we will get the header
        donors_snps = None
        # we will have a result per SNP
        results = []
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
                    # this could be made into a parallel process
                    #with Pool() as pool:
                    #    result = pool.map(self.filter_and_map_snp, [data_row, donor_offset, donors_snps, snp_id])
                    result = self.filter_and_map_snp(data_row, donor_offset, donors_snps, snp_id)
                    # add to the array
                    results.append(result)
                else:
                    warnings.warn(' '.join(['skipping', snp_id]))

            else:
                # need to set the header only once
                donors_snps = data_row[list(range(donor_offset, len(data_row), 1))]
                is_header = False
        # use nested comprehension to merge the list of lists
        results = [j for i in results for j in i]
        # turn the result into pandas dataframe
        results_pandas = pandas.DataFrame.from_dict(results)
        # do FDR correction
        results_pandas['BH'] = multitest.multipletests(pvals = results_pandas['p'].values.tolist(), method = 'fdr_bh')[1]
        # sort by corrected p
        results_pandas.sort_values(by = 'p', inplace = True)
        print('writing results to : ' + self.qtl_config.output_location)
        results_pandas.to_csv(self.qtl_config.output_location, sep = '\t', index = False)



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
