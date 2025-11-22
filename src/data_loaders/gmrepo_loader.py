import pandas as pd
import requests
import json
import time
from pandas import DataFrame

# COMMANDS URLS
BASE_URL = 'https://gmrepo.humangut.info/api/'
# get all the phenotypes
GET_ALL_PHENOTYPES = BASE_URL + 'get_all_phenotypes'
# get statistics for phenotype
GET_STATISTICS_BY_PROJECT_AND_PHENOTYPE = BASE_URL + 'getStatisticsByProjectsByMeshID'
GET_ASSOCIATED_SPECIES_WITH_PHENOTYPE = BASE_URL + 'getAssociatedSpeciesByMeshID'
GET_ASSOCIATED_GENERA_WITH_PHENOTYPE = BASE_URL + 'getAssociatedGeneraByMeshID'
GET_ASSOCIATED_PROJECTS_BY_PHENOTYPE = BASE_URL + 'getAssociatedProjectsByMeshID'
GET_ALL_PROJECTS = BASE_URL + 'getCuratedProjectsList'
GET_RUNS_COUNT_BY_PHENOTYPE = BASE_URL + 'countAssociatedRunsByPhenotypeMeshID'
GET_RUNS_BY_PHENOTYPE = BASE_URL + 'getAssociatedRunsByPhenotypeMeshIDLimit'
GET_MICROBIOME_ABUNDANCE_BY_PHENOTYPE_AND_TAXON = BASE_URL + 'getMicrobeAbundancesByPhenotypeMeshIDAndNCBITaxonID'
GET_ALL_MICROBES = BASE_URL + 'get_all_gut_microbes'
GET_PHENOTYPE_AND_ABUNDANCE_SUMMARY_BY_TAXON = BASE_URL + 'getPhenotypesAndAbundanceSummaryOfAAssociatedTaxon'
GET_ASSOCIATED_PHENOTYPES_AND_ABUNDANCE_OF_A_TAXON = BASE_URL + 'getAssociatedPhenotypesAndAbundancesOfATaxon'
GET_FULL_TAXONOMIC_PROFILE_BY_RUN_ID = BASE_URL + 'getFullTaxonomicProfileByRunID'
GET_MICROBE_ABUNDANCE_BY_PHENOTYPE_AND_PROJECT = BASE_URL + 'getMicrobeAbundancesByPhenotypeMeshIDAndProjectID'


class GMRepoLoader:

    @staticmethod
    def get_all_phenotypes() -> DataFrame:
        pheno_01 = requests.post(GET_ALL_PHENOTYPES, data={})
        pheno_01_cont = pheno_01.json().get('phenotypes')

        phenotypes = DataFrame(pheno_01_cont)
        return phenotypes

    @staticmethod
    def get_statistics_by_phenotype(phenotype_code: str) -> DataFrame:
        pheno_02_query = {'mesh_id': phenotype_code}
        pheno_02 = requests.post(GET_STATISTICS_BY_PROJECT_AND_PHENOTYPE, data=json.dumps(pheno_02_query))

        phenotype_stats = DataFrame(pheno_02.json())
        return phenotype_stats

    @staticmethod
    def get_species_associated_with_phenotype(phenotype_code: str) -> DataFrame:
        pheno_03_query = {'mesh_id': phenotype_code}
        pheno_03 = requests.post(GET_ASSOCIATED_SPECIES_WITH_PHENOTYPE, data=json.dumps(pheno_03_query))
        phenotyp_assoc_species = DataFrame(pheno_03.json())
        return phenotyp_assoc_species

    @staticmethod
    def get_genera_associated_with_phenotypes(phenotype_code: str) -> DataFrame:
        pheno_04_query = {'mesh_id': phenotype_code}
        pheno_04 = requests.post(GET_ASSOCIATED_GENERA_WITH_PHENOTYPE, data=json.dumps(pheno_04_query))
        phenotype_assoc_genera = DataFrame(pheno_04.json())
        return phenotype_assoc_genera

    @staticmethod
    def get_associated_projects_by_phenotypes(phenotype_code: str) -> DataFrame:
        pheno_05_query = {'mesh_id': phenotype_code}
        pheno_05 = requests.post(GET_ASSOCIATED_PROJECTS_BY_PHENOTYPE, data=json.dumps(pheno_05_query))
        phenotype_assoc_projects = DataFrame(pheno_05.json())

        return phenotype_assoc_projects

    @staticmethod
    def get_all_projects():
        query = {}
        content = requests.post(GET_ALL_PROJECTS, data=json.dumps(query))
        curated_pros = DataFrame(content.json())
        return curated_pros

    @staticmethod
    def get_runs_count_by_phenotype(phenotype_code: str) -> DataFrame:
        pheno_06_query = {'mesh_id': phenotype_code}
        pheno_06 = requests.post(GET_RUNS_COUNT_BY_PHENOTYPE, data=json.dumps(pheno_06_query))

        phenotype_nr_assoc_runs = DataFrame(pheno_06.json())
        runs_count = phenotype_nr_assoc_runs.loc[0, 'nr_assoc_runs']
        return runs_count

    @staticmethod
    def get_runs_by_phenotype(phenotype_code: str, skip: int = 0, limit: int = 100) -> DataFrame:
        pheno_07_query = {'mesh_id': phenotype_code, "skip": skip, "limit": limit}
        pheno_07 = requests.post(GET_RUNS_BY_PHENOTYPE, data=json.dumps(pheno_07_query))
        phenotype_a_page_of_assoc_runs = DataFrame(pheno_07.json())
        return phenotype_a_page_of_assoc_runs

    @staticmethod
    def get_microbiome_abundance_by_phenotype_and_taxon(phenotype_code: str, taxon_code: str) -> DataFrame:
        data_query = {'mesh_id': phenotype_code, "ncbi_taxon_id": taxon_code}
        data = requests.post(GET_MICROBIOME_ABUNDANCE_BY_PHENOTYPE_AND_TAXON, data=json.dumps(data_query))

        hist_data_for_phenotype = DataFrame(data.json().get('hist_data_for_phenotype'))
        return hist_data_for_phenotype

    @staticmethod
    def get_all_microbes_species() -> DataFrame:
        # get all species and genera that presented in >= 2 runs with median relative abundance >= 0.01%
        data = requests.post(GET_ALL_MICROBES, data={}).json()
        all_species = DataFrame(data.get('all_species'))
        return all_species

    @staticmethod
    def get_all_microbes_genera() -> DataFrame:
        # get all species and genera that presented in >= 2 runs with median relative abundance >= 0.01%
        data = requests.post(GET_ALL_MICROBES, data={}).json()
        all_genera = DataFrame(data.get('all_genus'))
        return all_genera

    @staticmethod
    def get_phenotype_and_abundance_summary_by_taxon(taxon_code: str) -> DataFrame:
        query = {"ncbi_taxon_id": taxon_code}
        data = DataFrame(requests.post(GET_PHENOTYPE_AND_ABUNDANCE_SUMMARY_BY_TAXON,
                                       data=json.dumps(query)).json().get('phenotypes_associated_with_taxon'))
        return data

    # This function does not work
    @staticmethod
    def get_associated_phenotypes_and_abundances_of_a_taxon(taxon_code: str) -> tuple[DataFrame, DataFrame, DataFrame]:
        query = {"ncbi_taxon_id": taxon_code}
        data = requests.post(GET_ASSOCIATED_PHENOTYPES_AND_ABUNDANCE_OF_A_TAXON, data=json.dumps(query)).json()

        phenotypes_associated_with_taxon = DataFrame(data.get("phenotypes_associated_with_taxon"))
        taxon = DataFrame(data.get("taxon"))
        density_data_groupped = DataFrame(data.get("density_data_groupped"))
        return phenotypes_associated_with_taxon, taxon, density_data_groupped

    @staticmethod
    def get_full_taxonomic_profile_by_run_id(run_id: str) -> tuple[dict, DataFrame, DataFrame]:
        query = {"run_id": run_id}
        data = requests.post(GET_FULL_TAXONOMIC_PROFILE_BY_RUN_ID, data=json.dumps(query)).json()
        run = data.get("run")
        species = DataFrame(data.get("species"))
        genus = DataFrame(data.get("genus"))
        return run, species, genus

    @staticmethod
    def get_microbe_abundances_by_phenotype_and_project(project_code: str,
                                                        phenotype_code: str = "") -> tuple[dict, dict, DataFrame]:
        query = {"project_id": project_code, "mesh_id": phenotype_code}

        # Query data
        data = requests.post(GET_MICROBE_ABUNDANCE_BY_PHENOTYPE_AND_PROJECT, data=json.dumps(query)).json()

        # Get project and disease information
        project = data.get("project_info")
        disease = data.get("disease_info")

        # Get abundance and meta data
        abundance_and_meta = DataFrame(data.get("abundance_and_meta_data"))
        return project, disease, abundance_and_meta

    def get_all_runs_by_phenotype(self, phenotype_code, inter_requests_delay: float = 0.5, step: int = 100) -> DataFrame:
        runs_count = self.get_runs_count_by_phenotype(phenotype_code)
        data = self.get_runs_by_phenotype(phenotype_code, 0, step)
        for i in range(step, runs_count, step):
            data = pd.concat([data, self.get_runs_by_phenotype(phenotype_code, i, step)])
            time.sleep(inter_requests_delay)
        return data

    def get_all_microbes_by_phenotype(self, phenotype_runs: DataFrame = None,
                                      phenotype_code: str = None, ) -> tuple[DataFrame, DataFrame]:
        if phenotype_runs is None:
            phenotype_runs = self.get_all_runs_by_phenotype(phenotype_code)
        data_genus = None
        data_species = None
        print('Runs:', len(phenotype_runs))
        count = 0
        for index, row in phenotype_runs.iterrows():
            run, species, genus = self.get_full_taxonomic_profile_by_run_id(row['run_id'])
            print(len(species))
            data_species = self.populate_dataframe_with_runs(data_species, species, run)
            data_genus = self.populate_dataframe_with_runs(data_genus, genus, run)
            count += 1
            print('Run:', count, 'out of', len(phenotype_runs))
            if data_species is not None:
                print("Species collected", len(data_species))
            if data_genus is not None:
                print("Genera collected:", len(data_genus))
        return data_genus, data_species

    def populate_dataframe_with_runs(self, all_data: DataFrame = None,
                                     data: DataFrame = None,
                                     run: dict = None) -> DataFrame:
        if not data.empty:
            data['run_id'] = run.get('run_id')
            data['project_id'] = run.get('project_id')
            data['host_age'] = run.get('host_age')
            data['sex'] = run.get('sex')
            data['BMI'] = run.get('BMI')
            data['country'] = run.get('country')
            if all_data is None:
                all_data = data.copy()
            else:
                all_data = pd.concat([all_data, data])
            return all_data
        else:
            return all_data


