
# external imports
import argparse

# internal exports
import QTLMapper

def parse_args():
    '''
    parse the command line arguments

    :return: a dictionary with keys and values like they would be created by argparse
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--snp_file_location', type = str, help = 'location of the snp file (string)')
    parser.add_argument('-p', '--probe_file_location', type = str, help = 'location of the probe file (string)')
    parser.add_argument('-c', '--covariates_file_location', type = str, help = 'location of the covariates file (string)')
    parser.add_argument('-sl', '--snp_positions_file_location', type = str, help = 'location of the snp positions file (string)')
    parser.add_argument('-pl', '--probe_positions_file_location', type = str, help = 'location of probe positions file (string)')
    parser.add_argument('-m', '--use_model', type = str, help = 'type of model to use (string)')
    parser.add_argument('-sc', '--confinements_snp_location', type = str, help = 'location snp confinement file (string)')
    parser.add_argument('-pc', '--confinements_probe_location', type = str, help = 'location of probe confinement file (string)')
    parser.add_argument('-cps', '--confinements_snp_probe_pairs_location', type = str, help = 'location snp-probe confinement file (string)')
    return parser


def get_test_args():
    test_args = {
        'snp_file_location' : '/Users/royoelen/hanze-master/2021/CeD_genotypes_adjusted27082018.txt',
        'probe_file_location' : '/Users/royoelen/hanze-master/2021/geuvadis_normalised_gene_expression_adjusted27082018.txt',
        'covariates_file_location' : None,
        'snp_positions_file_location' : '/Users/royoelen/hanze-master/2021/snp_locations_CeD_adjusted27082018.txt',
        'probe_positions_file_location' : '/Users/royoelen/hanze-master/2021/gene_locations.txt',
        'use_model' : 'linear',
        'confinements_snp_location' : None,
        'confinements_probe_location' : None,
        'confinements_snp_probe_pairs_location' : None,
        's' : '/Users/royoelen/hanze-master/2021/CeD_genotypes_adjusted27082018.txt',
        'p' : '/Users/royoelen/hanze-master/2021/geuvadis_normalised_gene_expression_adjusted27082018.txt',
        'c' : None,
        'sl' : None,
        'pl' : None,
        'm' : 'linear',
        'sc' : None,
        'pc' : None,
        'cps' : None
    }
    return test_args

@staticmethod
def args_to_QTLMapper(args):
    print(args)
    qtl_mapper = QTLMapper.QTLMapper(
        snp_file_location = args['snp_file_location'],
        probe_file_location = args['probe_file_location'],
        covariates_file_location = args['covariates_file_location'],
        snp_positions_file_location = args['snp_positions_file_location'],
        probe_positions_file_location = args['probe_positions_file_location'],
        use_model = args['use_model'],
        confinements_snp_location = args['confinements_snp_location'],
        confinements_probe_location = args['confinements_probe_location'],
        confinements_snp_probe_pairs_location = args['confinements_snp_probe_pairs_location']
    )
    return qtl_mapper

@staticmethod
def test():
    '''

    :return:
    '''
    # get the test arguments
    args = get_test_args()
    # create QTLMapper
    qtl_mapper = args_to_QTLMapper(args)
    # map
    qtl_mapper.perform_mapping()


# start
if __name__ == '__main__':
    '''
    application start
    '''
    test()
    #map()