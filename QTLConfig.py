

class QTLConfig:
    """
    Holder for configuration options in QTL mapping
    """

    def __init__(self,  snp_file_location=None, probe_file_location=None, output_location=None, covariates_file_location=None, covariates_to_use=None, snp_positions_file_location = None, probe_positions_file_location=None, use_model='linear', confinements_snp_location=None, confinements_probe_location=None, confinements_snp_probe_pairs_location=None, cis_dist=10000, maf=0.01):
        """
        constructor

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
        :param cis_dist: max distance from the probe start or stop, that a SNP is still considered within cis distance for
        :param maf: minimal minor allele frequency cutoff to consider a variant worth performing QTL mappinf for
        """
        self.snp_file_location = snp_file_location
        self.probe_file_location = probe_file_location
        self.output_location = output_location
        self.covariates_file_location = covariates_file_location
        self.covariates_to_use = covariates_to_use
        self.snp_positions_file_location = snp_positions_file_location
        self.probe_positions_file_location = probe_positions_file_location
        self.use_model = use_model
        self.confinements_snp_location = confinements_snp_location
        self.confinements_probe_location = confinements_probe_location
        self.confinements_snp_probe_pairs_location = confinements_snp_probe_pairs_location
        self.cis_dist = cis_dist
        self.maf = maf

