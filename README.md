# EuCanImage FHIR ETL Implementation

This repository contains the ETL implementation for EuCanImage, encouraging semantic interoperability of the clinical data obtained in the studies by transforming it into a machine-readable format following FHIR standards. This parser uses [FHIR Resources](https://github.com/nazrulworld/fhir.resources) in order to create the dictionaries following a FHIR compliant structure.
- Code Language is written in [Python 3.11](https://www.python.org/downloads/release/python-3110/).
- The outputs are JSON files compliant with [FHIR 4.3](https://hl7.org/fhir/R4B/) schemas.
- This script is specifically created for the Extract, Transform and Load implementation for EuCanImage, and will follow the structures obtained from the REDCap databases within the study. To create your own implementation in a different study, you may use the previously mentioned [FHIR Resources](https://github.com/nazrulworld/fhir.resources).

#### Data conversion process:
This code followed the structure to go through the following steps:
- Importing and transforming CSV with patient data
- Defining dictionaries for ontologies and functions to populate FHIR dictionaries
- Transforming dictionaries into FHIR resources
- Grouping FHIR resources into a defined bundle/envelope of resources
- Exporting as json file

#### Input & Output
- CSV file for each use case (CSV folder)
- JSON file following FHIR standards (OUTPUT folder)

## Installation and Guide
The first step is to clone or download the repository to your computer
```bash
git clone https://github.com/EGA-archive/EuCanImage-FHIR.git
```
#### Requirements
- Python 3.11.2
- [FHIR Resources](https://github.com/nazrulworld/fhir.resources) 6.5.0
- pandas 2.1.3
- numpy 1.26.2

In order to use these scripts, you will need to have access to [Python 3.11](https://www.python.org/downloads/release/python-3110/) in your systems.

To install the libraries used for this study, it can easily be done with `pip install`. The latest versions of each library should not cause any incompatibility.
```bash
pip install fhir.resources
pip install pandas
pip install numpy
```
### Instructions
The steps are the same on each Use Case, so we will be using Use Case 1 as an example for the steps to follow.

First of all, you will need to provide with a [CSV file](https://github.com/EGA-archive/EuCanImage-FHIR/blob/main/UC1_Hepatocellular_Carcinoma/CSV/UseCase1_testdata.csv) that follows the structure of the eCRF of the study. Each use case will have its own eCRF. Save the CSV file in the [CSV folder](https://github.com/EGA-archive/EuCanImage-FHIR/tree/main/UC1_Hepatocellular_Carcinoma/CSV) of the specific use case you will be using.

Next, in the beginning of each python file (For example, for Use Case 1 it would be [UC1-ETL.py](https://github.com/EGA-archive/EuCanImage-FHIR/blob/main/UC1_Hepatocellular_Carcinoma/UC1-ETL.py), you will need to change the variable `relative_path_csv` to change the name of the file matching the one of the input.
```bash
relative_path_csv = "/UC1_Hepatocellular_Carcinoma/CSV/UseCase1_testdata.csv"
```
Then, you can run the parser in the terminal, changing `<PATH-TO-FOLDER>` to the specific folder the parser is in, unless the terminal is run in the folder itself.
```bash
python <PATH-TO-FOLDER>/UC1-ETL.py
```
Once it is finished, you will have all of the parsed JSON files in the [OUTPUT](https://github.com/EGA-archive/EuCanImage-FHIR/tree/main/UC1_Hepatocellular_Carcinoma/OUTPUT) folder
