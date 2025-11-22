import os
import pickle

import pandas as pd
import pytest

from src.metabolome.metabolome_extractor import MetabolomeExtractor


@pytest.fixture(scope='session')  # one server to rule'em all
def metabolome_extractor():
    bacteria_files_map = pd.read_csv('../data/microbes_files_map_for_gmrepo2.csv', sep=';')
    metabolome_extractor = MetabolomeExtractor(bacteria_files_map)
    return metabolome_extractor


def test_get_metabolites(metabolome_extractor):
    pheno_data = pd.read_csv('../data/autistic_disorder_species.csv')
    metabolites = metabolome_extractor.compute_metabolites_from_gmrepo(pheno_data)
    print(metabolites)
    metabolites.to_csv('../data/metabolites_autistic_disorder_species.csv', index=False)


def test_create_metabolites_dictionary(metabolome_extractor):
    files_dict = {}
    files_list = os.listdir(metabolome_extractor.data_dir)
    index = 0
    for file in files_list:
        index += 1
        print(index, '/', len(files_list))
        files_dict[file] = metabolome_extractor.get_secreted_metabolites_from_file(file)
    with open('../AGORA-2.01/secreted_metabolites_by_file.pkl', 'wb') as fp:
        pickle.dump(files_dict, fp)


def test_check_if_genera_has_similar_metabolites(metabolome_extractor):
    with open(os.path.join('../AGORA-2.01/', 'secreted_metabolites_by_file.pkl'), 'rb') as fp:
        bacteria_metabolites_map = pickle.load(fp)
    sets = set()
    for bacteria_file in bacteria_metabolites_map:
        if "Haemophilus_pittmaniae_" in bacteria_file:
            sets.add(bacteria_metabolites_map[bacteria_file])

    print("\n\n\nNumber of species:", len(sets))
    print("\nTotal number of metabolites:", len(frozenset.union(*sets)))
    print("\nNumber of shared metabolites:", len(frozenset.intersection(*sets)))
    print("\nRatio between intersection and union", len(frozenset.intersection(*sets))/len(frozenset.union(*sets)))


