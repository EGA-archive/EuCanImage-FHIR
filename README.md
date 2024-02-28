Coding the use case JSON to make a template out of it.

#### Done in 5 steps:
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
