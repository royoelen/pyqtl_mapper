

class QTLConfig:


    def __init__(self,  snp_file_location=None, probe_file_location=None, covariates_file_location=None, snp_positions_file_location = None, probe_positions_file_location=None, use_model='linear', confinements_snp_location=None, confinements_probe_location=None, confinements_snp_probe_pairs_location=None, cis_dist=10000):
        '''

        :param snp_file_location:
        :param probe_file_location:
        :param covariates_file_location:
        :param snp_positions_file_location:
        :param probe_positions_file_location:
        :param use_model:
        :param confinements_snp_location:
        :param confinements_probe_location:
        :param confinements_snp_prob_pairs_location:
        '''
        self.snp_file_location = snp_file_location
        self.probe_file_location = probe_file_location
        self.covariates_file_location = covariates_file_location
        self.snp_positions_file_location = snp_positions_file_location
        self.probe_positions_file_location = probe_positions_file_location
        self.use_model = use_model
        self.confinements_snp_location = confinements_snp_location
        self.confinements_probe_location = confinements_probe_location
        self.confinements_snp_probe_pairs_location = confinements_snp_probe_pairs_location
        self.cis_dist = cis_dist

