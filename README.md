# IMPC: Individual Metabolite Production Capacity

This repository contains code for computing Individual Metabolite Production Capacity (IMPC), a framework designed for quantitative analysis of metabolite biosynthetic potential at the individual level. For more details please refer to the paper 
*"Machine Learning model to identify gut microbiome-derived metabolites as potential biomarkers Autism Spectrum Disorder: a pilot study"*

## About

IMPC provides a set of computational tools and workflows tailored for the estimation of metabolite production capacity based on genomic, metagenomic, or metabolomic data. This toolkit is applicable to studies in microbiome research, personalized medicine, and systems biology.

## Features

- Modular Python and Jupyter-based workflows.
- Methods for preprocessing, analysis, and visualization of IMPC metrics.
- Interface for customizable pipelines using your own datasets.
- Open-source: freely reusable under the MIT License.[1]

## Getting Started

### Prerequisites

- Python 3.x
- Jupyter Notebook
- Required Python packages (see your `requirements.txt` or specify here if available).

### Installation

Clone the repository:
```bash
git clone https://github.com/dec0der0b/IMPC.git
cd IMPC
```

Install dependencies:
```bash
pip install -r requirements.txt
```
or install required packages manually if no requirements file is provided.
Extract the AGORA-2.01.zip into the main folder.
### Usage

- Open the main notebook in Jupyter:
  ```bash
  jupyter notebook
  ```
- Follow the code cells and documentation in the notebook to input your own data and compute IMPC metrics.

## Code Structure

- **Jupyter Notebooks:** Main scripts for running the IMPC computation. It includes analyses for all ages, under 25 and above 25
- **Python Files:** Supporting functions and modules.
- **Data:** (Optional) Place for input datasets and example results.

## License

This project is licensed under the MIT License.

## Contributing

Please submit issues or pull requests for improvements, bugfixes, or new feature suggestions.

***
 
