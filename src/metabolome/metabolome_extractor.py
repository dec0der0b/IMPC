import pandas as pd
import cobra
import os
import glob
from cobra.core import Model
import pickle


from pandas import DataFrame

MICROBIOME_DATA_DIR = '../../AGORA-2.01'


class MetabolomeExtractor:
    def __init__(self, bacteria_files_map: DataFrame, data_dir: str = MICROBIOME_DATA_DIR):
        self.data_dir = os.path.join(data_dir, 'reconstructions', 'mat')
        self.bacteria_files_map = bacteria_files_map
        self.bacteria_metabolites_map = None
        try:
            with open(os.path.join(data_dir, 'secreted_metabolites_by_file.pkl'), 'rb') as fp:
                self.bacteria_metabolites_map = pickle.load(fp)
        except:
            print('Working without bacteria_metabolites_map')

    def compute_metabolites_extensive(self, species: DataFrame):
        metabolites_dict = dict()
        for index, row in species.iterrows():
            bacter_metabolites = set()
            bacter_id = row['ncbi_taxon_id']
            try:
                files = self.bacteria_files_map[self.bacteria_files_map['ncbi_taxon_id'] == bacter_id].iloc[0][
                        'files'].strip('][').split(', ')
            except:
                print('bacter:', row['scientific_name'], 'is not in bacteria list.')
                continue
            for file in files:
                if file == '':
                    continue
                file = file.replace("'", "").split('\\')[-1]
                metabolites = self.get_secreted_metabolites_from_file(file)
                bacter_metabolites.update(metabolites)
            for metabolite in bacter_metabolites:
                if metabolite in metabolites_dict.keys():
                    metabolites_dict[metabolite][0] += row['relative_abundance']
                    metabolites_dict[metabolite][1] += '/' + row['scientific_name']
                else:
                    metabolites_dict[metabolite] = [row['relative_abundance'], row['scientific_name']]
        if not bool(metabolites_dict):
            return None
        metabolites_df = pd.DataFrame.from_dict(metabolites_dict, orient='index')
        metabolites_df = metabolites_df.reset_index()
        metabolites_df.columns = ['metabolite', 'relative_abundance', 'bacteria']
        return metabolites_df

    def compute_metabolites_restrictive(self, species: DataFrame):
        metabolites_dict = dict()
        for index, row in species.iterrows():
            bacter_metabolites = set()
            bacter_id = row['ncbi_taxon_id']
            try:
                files = self.bacteria_files_map[self.bacteria_files_map['ncbi_taxon_id'] == bacter_id].iloc[0][
                    'files'].strip('][').split(', ')
            except:
                print('bacter:', row['scientific_name'], 'is not in bacteria list.')
                continue
            for file in files:
                if file == '':
                    continue
                file = file.replace("'", "").split('\\')[-1]
                metabolites = self.get_secreted_metabolites_from_file(file)
                if len(bacter_metabolites) == 0:
                    bacter_metabolites.update(metabolites)
                else:
                    bacter_metabolites = bacter_metabolites.intersection(metabolites)
            for metabolite in bacter_metabolites:
                if metabolite in metabolites_dict.keys():
                    metabolites_dict[metabolite][0] += row['relative_abundance']
                    metabolites_dict[metabolite][1] += '/' + row['scientific_name']
                else:
                    metabolites_dict[metabolite] = [row['relative_abundance'], row['scientific_name']]
        if not bool(metabolites_dict):
            return None
        metabolites_df = pd.DataFrame.from_dict(metabolites_dict, orient='index')
        metabolites_df = metabolites_df.reset_index()
        metabolites_df.columns = ['metabolite', 'relative_abundance', 'bacteria']
        return metabolites_df

    def compute_metabolites_from_gmrepo(self, df) -> DataFrame:
        metabolites_dict = dict()
        for index, row in df.iterrows():
            bacter_metabolites = set()
            bacter_id = row['ncbi_taxon_id']
            try:
                files = self.bacteria_files_map[self.bacteria_files_map['ncbi_taxon_id'] == bacter_id].iloc[0]['files'].strip('][').split(', ')
            except:
                print('bacter:', row['scientific_name'], 'is not in bacteria list.')
                continue
            for file in files:
                if file == '':
                    continue
                file = file.replace("'", "").split('\\')[-1]
                metabolites = self.get_secreted_metabolites_from_file(file)
                bacter_metabolites.update(metabolites)
            for metabolite in bacter_metabolites:
                if metabolite in metabolites_dict.keys():
                    metabolites_dict[metabolite][0] += row['abus_mean']
                    metabolites_dict[metabolite][1] += row['abus_median']
                    metabolites_dict[metabolite][2] += row['abus_mean'] - row['abus_sd']
                    metabolites_dict[metabolite][3] += row['abus_mean'] + row['abus_sd']
                    metabolites_dict[metabolite][4] += '/' + row['scientific_name']
                else:
                    metabolites_dict[metabolite] = [row['abus_mean'], row['abus_median'],
                                                    row['abus_mean'] - row['abus_sd'],
                                                    row['abus_mean'] + row['abus_sd'],
                                                    row['scientific_name']]
        metabolites_df = pd.DataFrame.from_dict(metabolites_dict, orient='index')
        metabolites_df = metabolites_df.reset_index()
        metabolites_df.columns = ['metabolite', 'abus_mean', 'abus_median', 'abus_min', 'abus_max', 'bacteria']
        return metabolites_df

    @staticmethod
    def get_secreted_metabolites_from_model(model: Model) -> set:
        try:
            secretions = model.summary().secretion_flux
            result = secretions['metabolite'].tolist()
            return set(result)
        except:
            print('error on ' + model.name)
            return set()

    @staticmethod
    def get_uptake_metabolites_from_model(model: Model) -> set:
        try:
            secretions = model.summary().uptake_flux
            result = secretions['metabolite'].tolist()
            return set(result)
        except:
            print('error on ' + model.name)
            return set()

    @staticmethod
    def replace_metabolites_id_with_name(secreted_metabolites, metabolites):
        result = set()
        for metabolite in metabolites:
            if metabolite.id in secreted_metabolites:
                result.add(metabolite.name)
        return frozenset(result)

    def get_secreted_metabolites_from_file(self, filename: str) -> frozenset:
        if self.bacteria_metabolites_map is not None:
            return self.bacteria_metabolites_map[filename]
        model = cobra.io.load_matlab_model(os.path.join(self.data_dir, filename))
        secreted_metabolites = self.get_secreted_metabolites_from_model(model)
        metabolites = model.metabolites
        final_metabolites = self.replace_metabolites_id_with_name(secreted_metabolites, metabolites)
        return final_metabolites

    def get_uptake_metabolites_from_file(self, filename: str) -> frozenset:
        model = cobra.io.read_sbml_model(os.path.join(self.data_dir, filename.replace('.xml', '.mat')))
        secreted_metabolites = self.get_uptake_metabolites_from_model(model)
        metabolites = model.metabolites
        final_metabolites = self.replace_metabolites_id_with_name(secreted_metabolites, metabolites)
        return final_metabolites

    def add_filename_to_bacter(self, row) -> list:
        bacter_name_parts = row['scientific_name'].split()
        if len(bacter_name_parts) == 2:
            bacter_name = bacter_name_parts[0].replace('[', '').replace(']', '') + '_' + bacter_name_parts[1]
        filenamesList = glob.glob(os.path.join(self.data_dir, bacter_name + '*.xml'))
        row['files'] = filenamesList
        return row
