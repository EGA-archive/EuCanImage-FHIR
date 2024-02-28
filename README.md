ETL implementation for EuCanImage, encouraging semantic interoperability of the clinical data obtained in the studies by transforming it into a machine-readable format following FHIR standards.

#### 5 steps:
- Importing and transforming CSV with patient data
- Defining dictionaries for ontologies and functions to populate FHIR dictionaries
- Transforming dictionaries into FHIR resources
- Grouping FHIR resources into a defined bundle/envelope of resources
- Exporting as json file

### Requirements:
- Python 3.11.2
- [FHIR Resources](https://github.com/nazrulworld/fhir.resources) 6.5.0
- pandas 2.1.3
- numpy 1.26.2

#### Input & Output
- CSV file for each use case (CSV folder)
- JSON file following FHIR standards (OUTPUT folder)