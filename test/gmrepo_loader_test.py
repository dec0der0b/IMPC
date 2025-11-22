import os
import glob

import pytest
import pandas as pd

from src.data_loaders.gmrepo_loader import GMRepoLoader

ALL_PHENOTYPES_DICT = {
    'Colitis Ulcerative': 'D003093',
    'Parkinson Disease': 'D010300',
    'Irritable Bowel Syndrome': 'D043183',
    'Diarrhea': 'D003967',
    'Cystic Fibrosis': 'D003550',
    'Constipation': 'D003248',
    'Migraine Disorders': 'D008881',
    'Autoimmune Diseases': 'D001327',
    'Inflammatory Bowel Diseases': 'D015212',
    'Bipolar Disorder': 'D001714',
    'Intestinal Diseases': 'D007410',
    'Depression': 'D003863',
    'Schizophrenia': 'D012559',
    'Attention Deficit Disorder with Hyperactivity': 'D001289',
    'Arthritis, Rheumatoid': 'D001172',
    'Alzheimer Disease': 'D000544',
    'Celiac Disease': 'D002446',
    'Gastroesophageal Reflux': 'D005764',
    'Chronic Fatigue Syndrome': 'D015673',
    'Metabolic Syndrome': 'D024821',
    'Cognitive Dysfunction': 'D060825',
    'Rett Syndrome': 'D015518',
    'Health': 'D006262',
    'Autism Spectrum Disorder': 'D000067877',
    'Autistic Disorder': 'D001321',
    'Epilepsy': 'D004827'
}


@pytest.fixture(scope='session')  # one server to rule'em all
def gmrepo_loader():
    gmrepo_loader = GMRepoLoader()
    return gmrepo_loader


def test_get_all_phenotypes(gmrepo_loader):
    phenotypes = gmrepo_loader.get_all_phenotypes()
    print(phenotypes.head())
    assert len(phenotypes) == 132
    phenotypes.to_csv('../data/phenotypes.csv', index=False)


def test_get_statistics_by_phenotype(gmrepo_loader):
    statistics = gmrepo_loader.get_statistics_by_phenotype('D006262')
    print(statistics)
    assert statistics.loc['term', 'metadata'] == 'Health'


def test_get_species_associated_with_phenotype(gmrepo_loader):
    species = gmrepo_loader.get_species_associated_with_phenotype('D006262')
    print(species.head())
    assert species.loc[0, 'disease'] == 'D006262'


def test_get_genera_associated_with_phenotypes(gmrepo_loader):
    genera = gmrepo_loader.get_species_associated_with_phenotype('D006262')
    print(genera.head())
    assert genera.loc[0, 'disease'] == 'D006262'


def test_get_associated_projects_by_phenotypes(gmrepo_loader):
    genera = gmrepo_loader.get_associated_projects_by_phenotypes('D006262')
    print(genera.head())
    assert genera.loc[0, 'mesh_id'] == 'D006262'


def test_get_all_projects(gmrepo_loader):
    projects = gmrepo_loader.get_all_projects()
    print(projects.head())
    assert projects.loc[0, 'project_id'] == 'PRJDB4176'


def test_get_runs_count_by_phenotype(gmrepo_loader):
    runs_count = gmrepo_loader.get_runs_count_by_phenotype('D006262')
    print(runs_count)
    assert runs_count == 34019


def test_get_runs_by_phenotype(gmrepo_loader):
    runs = gmrepo_loader.get_runs_by_phenotype('D006262', 0, 10)
    print(runs.head())
    assert runs.loc[0, 'project_id'] == 'PRJDB3601'


def test_get_microbiome_abundance_by_phenotype_and_taxon(gmrepo_loader):
    microbiome_abundance = gmrepo_loader.get_microbiome_abundance_by_phenotype_and_taxon('D003093', '40520')
    print(microbiome_abundance.head())
    assert round(microbiome_abundance.loc[0, 'cumpct'], 2) == 89.78


def test_get_all_microbes(gmrepo_loader):
    microbes = gmrepo_loader.get_all_microbes()
    print(microbes.head())
    assert microbes.loc[0, 'name'] == 'Azorhizobium caulinodans'


def test_get_all_microbes_genera(gmrepo_loader):
    microbes = gmrepo_loader.get_all_microbes_genera()
    print(microbes.head())
    assert microbes.loc[0, 'name'] == 'Azorhizobium'


def test_get_phenotype_and_abundance_summary_by_taxon(gmrepo_loader):
    phen_and_abundance = gmrepo_loader.get_phenotype_and_abundance_summary_by_taxon('40520')
    print(phen_and_abundance.head())
    assert phen_and_abundance.loc[0, 'disease'] == 'D025241'


def test_get_associated_phenotypes_and_abundances_of_a_taxon(gmrepo_loader):
    phenotypes_associated_with_taxon, taxon, density_data_groupped = gmrepo_loader.\
        get_associated_phenotypes_and_abundances_of_a_taxon('40520')
    print(phenotypes_associated_with_taxon.head())
    print(taxon)
    print(density_data_groupped)
    assert density_data_groupped.loc[0, 'disease'] == 'D025241'


def test_get_full_taxonomic_profile_by_run_id(gmrepo_loader):
    run, species, genus = gmrepo_loader.get_full_taxonomic_profile_by_run_id("SRR7811020")
    print(run)
    print(species.head())
    print(genus.head())
    assert genus.loc[0, 'scientific_name'] == 'Bacteroides'


def test_get_microbe_abundances_by_phenotype_and_project(gmrepo_loader):
    project, disease, abundance_and_meta = gmrepo_loader.\
        get_microbe_abundances_by_phenotype_and_project("PRJNA489760", "D006262")
    print(project)
    print(disease)
    print(abundance_and_meta.head())
    assert isinstance(disease, dict)
    assert abundance_and_meta.loc[0, 'scientific_name'] == 'Bifidobacterium'

    project, disease, abundance_and_meta = gmrepo_loader. \
        get_microbe_abundances_by_phenotype_and_project("PRJNA489760", "")
    print(project)
    print(disease)
    print(abundance_and_meta.head())
    assert isinstance(disease, list)
    assert abundance_and_meta.loc[0, 'scientific_name'] == 'Bifidobacterium'


def test_get_all_microbes_by_phenotype(gmrepo_loader):
    data_genus, data_species = gmrepo_loader.get_all_microbes_by_phenotype(phenotype_code='D001321')
    if data_genus is not None:
        data_genus.to_csv('../data/Autistic_Disorder_genus.csv', index=False)
    if data_species is not None:
        data_species.to_csv('../data/Autistic_Disorder_species.csv', index=False)


def test_find_files_for_microbe(gmrepo_loader):
    microbes = gmrepo_loader.get_all_microbes_species()
    data_dir = os.path.join('../AGORA-2.01', 'reconstructions', 'mat')

    def add_filename_to_bacter(row) -> list:
        bacter_name_parts = row['name'].split()
        bacter_name = bacter_name_parts[0].replace('[', '').replace(']', '') + '_' + bacter_name_parts[1].replace('.', '_')
        filenamesList = glob.glob(os.path.join(data_dir, bacter_name + '*.mat'))
        row['files'] = filenamesList
        return row

    microbes = microbes.apply(lambda row: add_filename_to_bacter(row), axis=1)
    microbes.to_csv('../data/microbes_species_agora2.csv', index=False)


def test_find_files_for_microbe_genera(gmrepo_loader):
    microbes = gmrepo_loader.get_all_microbes_genera()
    data_dir = os.path.join('../AGORA-2.01', 'reconstructions', 'mat')

    def add_filename_to_bacter(row) -> list:
        bacter_name_parts = row['name'].split()
        bacter_name = bacter_name_parts[0].replace('[', '').replace(']', '') + '_'
        filenamesList = glob.glob(os.path.join(data_dir, bacter_name + '*.mat'))
        row['files'] = [filename.split('\\')[-1] for filename in filenamesList]
        return row

    microbes = microbes.apply(lambda row: add_filename_to_bacter(row), axis=1)
    microbes.to_csv('../data/microbes_genera_agora2.csv', index=False)


def test_get_all_microbes_by_all_phenotype(gmrepo_loader):
    for phenotype, code in ALL_PHENOTYPES_DICT.items():
        print(phenotype, code)
        data_genus, data_species = gmrepo_loader.get_all_microbes_by_phenotype(phenotype_code=code)
        if data_genus is not None:
            data_genus.to_csv(f'../data/all_phenotypes/{phenotype}_genus.csv', index=False)
        if data_species is not None:
            data_species.to_csv(f'../data/all_phenotypes/{phenotype}_species.csv', index=False)









