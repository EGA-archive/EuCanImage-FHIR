# %%
#Importing modules
import os
import pandas as pd
import numpy as np
import math
import uuid
from fhir.resources.bundle import Bundle
from fhir.resources.patient import Patient
from fhir.resources.condition import Condition
from fhir.resources.observation import Observation
from fhir.resources.medicationadministration import MedicationAdministration
import json

# %%
#Import data from CSV. Change formats to readable ones by dictionaries.
#absolute_path = os.path.dirname(os.getcwd())
absolute_path = os.getcwd()
relative_path_csv = "/UC6&8/CSV/UseCase68_testdata.csv"
csvpath = absolute_path+relative_path_csv
data_csv = pd.read_csv(csvpath)
data_csv = data_csv.astype({
    "sex":"Int64",
    "patientclass":"Int64",
    "menop":"Int64",
    "n_preg":"Int64",
    "lactation":"Int64",
    "symptoms":"Int64",
    "screening":"Int64",
    "famhisto_b":"Int64",
    "famhisto_o":"Int64",
    "laterality":"Int64",
    "bcstproxy":"Int64",
    "bcst_pam50":"Int64",
    "histopat":"Int64",
    "histopat_2":"Int64",
    "histopat_3":"Int64",
    "tnm_ct":"Int64",
    "tnm_cn":"Int64",
    "tnm_cm":"Int64",
    "grade":"Int64",
    "dcis":"Int64",
    "her2ihc":"Int64",
    "her2fish":"Int64",
    "laterality_2":"Int64",
    "bcstproxy_2":"Int64",
    "bcst_pam50_2":"Int64",
    "histopat_1_2":"Int64",
    "histopat_2_2":"Int64",
    "histopat_3_2":"Int64",
    "tnm_ct_2":"Int64",
    "tnm_cn_2":"Int64",
    "tnm_cm_2":"Int64",
    "grade_2":"Int64",
    "dcis_2":"Int64",
    "her2ihc_2":"Int64",
    "her2fish_2":"Int64",
    "hcontrac":"Int64",
    "hormtherapy":"Int64",
    "brca1":"Int64",
    "brca2":"Int64",
    "palb2":"Int64",
    "chek2":"Int64" 
})

#Convert table to dictionary containing all variables as lists
d_csv = {}
for i in range(len(data_csv.columns)):
    d_csv[data_csv.columns[i]+"_csv"] = data_csv[data_csv.columns[i]].tolist()

# %%
sex_code= {"system": "http://loinc.org","code": "76689-9","display": "Sex assigned at birth"}
sex_dict = {0: {"system":"http://loinc.org","code":"LA2-8","display":"Male"},
            1: {"system":"http://loinc.org","code":"LA3-6","display":"Female"}}
patientclass_dict = {0:{"system":"http://snomed.info/sct","code":"93796005","display":"Primary malignant neoplasm of female breast"},
                     1:{"system":"http://snomed.info/sct","code":"3898006","display":"Neoplasm, benign (morphologic abnormality)"},
                     2:{"system":"http://snomed.info/sct","code":"17621005","display":"Normal (qualifier value)"}}
                    #2:{"system":"http://snomed.info/sct","code":"168749009","display":"Mammography normal (finding)"}}
bodysite_code = {"system":"http://snomed.info/sct","code":"76752008","display":"Breast structure (body structure)"}
menop_code = {"system":"http://snomed.info/sct","code":"161712005","display":"Menopause, function (observable entity)"}
menop_dict = {0:{"system":"http://snomed.info/sct","code":"289904000","display":"Menopause absent (finding)"},
              1:{"system":"http://snomed.info/sct","code":"289903006","display":"Menopause present (finding)"}}
n_preg_code = {"system":"http://snomed.info/sct","code":"161732006","display":"Gravida (observable entity)"}
lactation_code = {"system":"http://snomed.info/sct","code":"169741004","display":"Breast fed (finding)"}
symptoms_code = {"system":"http://snomed.info/sct","code":"198116001","display":"Breast signs and symptoms (finding)"}
screening_code = {"system":"http://snomed.info/sct","code":"609223006","display":"Magnetic resonance imaging of breast for screening for malignant neoplasm (procedure)"}
famhisto_b_code = {"system":"http://snomed.info/sct","code":"416471007:64572001=372064008","display":"Family history of clinical finding where Disease = Malignant neoplasm of female breast"}
famhisto_o_code = {"system":"http://snomed.info/sct","code":"416471007:64572001=363443007","display":"Family history of clinical finding where Disease = Malignant tumour of ovary"}
famhisto_dict = {0:{"system":"http://snomed.info/sct","code":"260413007","display":"None (qualifier value)"},
                 1:{"system":"http://snomed.info/sct","code":"264500008","display":"First degree (qualifier value)"},
                 2:{"system":"http://snomed.info/sct","code":"263868004","display":"Second degree (qualifier value)"},
                 3:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
laterality_dict = {0:{"system":"http://snomed.info/sct","code":"80248007","display":"Left breast structure (body structure)"},
                   1:{"system":"http://snomed.info/sct","code":"73056007","display":"Right breast structure (body structure)"}}
bcstproxy_code = {"system":"http://snomed.info/sct","code":"372064008+260837004","display":"Malignant neoplasm of female breast, and Subtype"}
bcst_pam50_code = {"system":"http://snomed.info/sct","code":"372064008+260837004","display":"Malignant neoplasm of female breast, and Subtype"}
bc_dict={0:{"system":"http://snomed.info/sct","code":"706970001","display":"Triple negative malignant neoplasm of breast (disorder)"},
         1:[{"code":{"coding":[{"system":"http://loinc.org","code":"16112-5","display":"Estrogen receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"416053008","display":"Estrogen receptor positive tumor (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"16113-3","display":"Progesterone receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"416561008","display":"Progesterone receptor positive tumor (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"48676-1","display":"HER2 [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"431396003","display":"Human epidermal growth factor 2 negative carcinoma of breast (disorder)"}]}},
            {"code":{"coding":[{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C123557","display":"Ki67 Measurement"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C162076","display":"MKI67 Negative"}]}}],
         2:[{"code":{"coding":[{"system":"http://loinc.org","code":"16112-5","display":"Estrogen receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"416053008","display":"Estrogen receptor positive tumor (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"16113-3","display":"Progesterone receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"416561008","display":"Progesterone receptor positive tumor (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"48676-1","display":"HER2 [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"431396003","display":"Human epidermal growth factor 2 negative carcinoma of breast (disorder)"}]}},
            {"code":{"coding":[{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C123557","display":"Ki67 Measurement"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C146686","display":"MKI67 Positive"}]}}],
         3:[{"code":{"coding":[{"system":"http://loinc.org","code":"16112-5","display":"Estrogen receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"416053008","display":"Estrogen receptor positive tumor (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"16113-3","display":"Progesterone receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"416561008","display":"Progesterone receptor positive tumor (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"48676-1","display":"HER2 [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"427685000","display":"Human epidermal growth factor 2 positive carcinoma of breast (disorder)"}]}},
            {"code":{"coding":[{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C123557","display":"Ki67 Measurement"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C146686","display":"MKI67 Positive"}]}}],
         4:[{"code":{"coding":[{"system":"http://loinc.org","code":"16112-5","display":"Estrogen receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"441117001","display":"Estrogen receptor negative neoplasm (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"16113-3","display":"Progesterone receptor [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"441118006","display":"Progesterone receptor negative neoplasm (disorder)"}]}},
            {"code":{"coding":[{"system":"http://loinc.org","code":"48676-1","display":"HER2 [Interpretation] in Tissue"}]},
             "valueCodeableConcept":{"coding":[{"system":"http://snomed.info/sct","code":"427685000","display":"Human epidermal growth factor 2 positive carcinoma of breast (disorder)"}]}}]}
bc_note_dict = {0:"Triple negative",
                1:"Luminal A",
                2:"Luminal B- HER2 negative",
                3:"HER2 positive (enriched)/HR positive (luminal)",
                4:"HER2 positive(enriched)/HR negative (non-luminal)"}
histopat_dict = {0:{"system":"http://snomed.info/sct","code":"3898006","display":"Neoplasm, benign (morphologic abnormality)"},
                 1:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8500/3","display":"DUCT CARCINOMA/ Invasive carcinoma of no special type"},
                 2:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8522/3","display":"LOBULAR AND OTHER DUCTAL CA. / Infiltrating duct and lobular carcinoma"},
                 3:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8523/3","display":"LOBULAR AND OTHER DUCTAL CA. / Infiltr. duct mixed with other types of carcinoma"},
                 4:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8520/3","display":"LOBULAR AND OTHER DUCTAL CA. / Lobular carcinoma, NOS"},
                 5:{"system":"http://snomed.info/sct","code":"713609000:410657003=74964007","display":"Invasive carcinoma of breast where Type = Other"},
                 6:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8500/2","display":"DUCT CARCINOMA/ Intraductal carcinoma, noninfiltrating, NOS"},
                 7:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8520/2","display":"LOBULAR AND OTHER DUCTAL CA. / Lobular carcinoma insitu"},
                 8:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8522/2","display":"LOBULAR AND OTHER DUCTAL CA. / Intraductal and lobular in situ carcinoma"},
                 9:{"system":"http://snomed.info/sct","code":"1187138006","display":"Carcinoma in situ (morphologic abnormality)"}}
tnm_ct_code = {"system":"http://snomed.info/sct","code":"1222585009","display":"American Joint Committee on Cancer clinical T category allowable value (qualifier value)"}
tnm_ct_dict = {0:{"system":"http://snomed.info/sct","code":"1222604002","display":"cTX"},
               1:{"system":"http://snomed.info/sct","code":"1228882005","display":"cT0"},
               2:{"system":"http://snomed.info/sct","code":"1228884006","display":"cTis"},
               3:{"system":"http://snomed.info/sct","code":"1228885007","display":"cTis(DCIS)"},
               4:{"system":"http://snomed.info/sct","code":"1228888009","display":"cTis(Paget)"},
               5:{"system":"http://snomed.info/sct","code":"1228889001","display":"cT1"},
               6:{"system":"http://snomed.info/sct","code":"1228892002","display":"cT1a"},
               7:{"system":"http://snomed.info/sct","code":"1228895000","display":"cT1b"},
               8:{"system":"http://snomed.info/sct","code":"1228899006","display":"cT1c"},
               9:{"system":"http://snomed.info/sct","code":"1228891009","display":"cT1mi"},
               10:{"system":"http://snomed.info/sct","code":"1228929004","display":"cT2"},
               11:{"system":"http://snomed.info/sct","code":"1228938002","display":"cT3"},
               12:{"system":"http://snomed.info/sct","code":"1228944003","display":"cT4"},
               13:{"system":"http://snomed.info/sct","code":"1228945002","display":"cT4a"},
               14:{"system":"http://snomed.info/sct","code":"1228946001","display":"cT4b"},
               15:{"system":"http://snomed.info/sct","code":"1228947005","display":"cT4c"},
               16:{"system":"http://snomed.info/sct","code":"1228948000","display":"cT4d"}}
tnm_cn_code = {"system":"http://snomed.info/sct","code":"1222588006","display":"American Joint Committee on Cancer clinical N category allowable value (qualifier value)"}
tnm_cn_dict = {0:{"system":"http://snomed.info/sct","code":"1229966003","display":"cNX"},
               1:{"system":"http://snomed.info/sct","code":"1229967007","display":"cN0"},
               2:{"system":"http://snomed.info/sct","code":"1229973008","display":"cN1"},
               3:{"system":"http://snomed.info/sct","code":"1229978004","display":"cN2"},
               4:{"system":"http://snomed.info/sct","code":"1229981009","display":"cN2a"},
               5:{"system":"http://snomed.info/sct","code":"1229982002","display":"cN2b"},
               6:{"system":"http://snomed.info/sct","code":"1229984001","display":"cN3"},
               7:{"system":"http://snomed.info/sct","code":"1229985000","display":"cN3a"},
               8:{"system":"http://snomed.info/sct","code":"1229986004","display":"cN3b"},
               9:{"system":"http://snomed.info/sct","code":"1229987008","display":"cN3c"}}
tnm_cm_code = {"system":"http://snomed.info/sct","code":"1222587001","display":"American Joint Committee on Cancer clinical M category allowable value (qualifier value)"}
tnm_cm_dict = {0:{"system":"http://snomed.info/sct","code":"1229901006","display":"cM0"},
               1:{"system":"http://snomed.info/sct","code":"1229903009","display":"cM1"},
               2:{"system":"http://snomed.info/sct","code":"1229904003","display":"cM1a"},
               3:{"system":"http://snomed.info/sct","code":"1229907005","display":"cM1b"},
               4:{"system":"http://snomed.info/sct","code":"1229910003","display":"cM1c"},
               5:{"system":"http://snomed.info/sct","code":"1229913001","display":"cM1d"}}
grade_code = {"system":"http://snomed.info/sct","code":"258244004","display":"Tumor histopathological grade status values (tumor staging)"}
grade_dict = {0:{"system":"http://snomed.info/sct","code":"1228845001","display":"American Joint Committee on Cancer grade GX (qualifier value)"},
              1:{"system":"http://snomed.info/sct","code":"1228848004","display":"American Joint Committee on Cancer grade G1 (qualifier value)"},
              2:{"system":"http://snomed.info/sct","code":"1228850007","display":"American Joint Committee on Cancer grade G2 (qualifier value)"},
              3:{"system":"http://snomed.info/sct","code":"1228851006","display":"American Joint Committee on Cancer grade G3 (qualifier value)"}}
dcis_code = {"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C116341","display":"DCIS Score"}
dcis_dict = {0:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C9457","display":"DCIS Grade 1"},
             1:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C9456","display":"DCIS Grade 2"},
             2:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C7949","display":"DCIS Grade 3"}}
er_code = {"system":"http://loinc.org","code":"85329-1","display":"Cells.estrogen receptor/100 cells in Breast cancer specimen by Immune stain"}
pr_code = {"system":"http://loinc.org","code":"85325-9","display":"Cells.progesterone receptor/100 cells in Breast cancer specimen by Immune stain"}
her2ihc_code = {"system":"http://loinc.org","code":"85319-2","display":"HER2 [Presence] in Breast cancer specimen by Immune stain"}
her2ihc_dict = {0:{"system":"http://loinc.org","code":"LA6111-4","display":"0"},
                1:{"system":"http://loinc.org","code":"LA11841-6","display":"1+"},
                2:{"system":"http://loinc.org","code":"LA11842-4","display":"2+"},
                3:{"system":"http://loinc.org","code":"LA11843-2","display":"3+"}}
ki67_code = {"system":"http://loinc.org","code":"85330-9","display":"Cells.Ki-67 nuclear Ag/100 cells in Breast cancer specimen by Immune stain"}
her2fish_code = {"system":"http://loinc.org","code":"85318-4","display":"ERBB2 gene duplication [Presence] in Breast cancer specimen by FISH"}
hcontrac_code = {"system":"http://snomed.info/sct","code":"1237404009","display":"Uses hormone method of contraception (finding)"}
hormtherapy_code = {"system":"http://snomed.info/sct","code":"266717002","display":"Hormone replacement therapy (procedure)"}
brca1_code = {"system":"http://snomed.info/sct","code":"405823003","display":"BRCA1 mutation carrier detection test (procedure)"}
brca1_dict = {0:{"system":"http://snomed.info/sct","code":"412734009","display":"BRCA1 gene mutation positive (finding)"},
              1:{"system":"http://snomed.info/sct","code":"412736006","display":"BRCA1 gene mutation negative (finding)"},
              2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown (qualifier value)"}}
brca2_code = {"system":"http://snomed.info/sct","code":"405826006","display":"BRCA2 mutation carrier detection test (procedure)"}
brca2_dict = {0:{"system":"http://snomed.info/sct","code":"412738007","display":"BRCA2 gene mutation positive (finding)"},
              1:{"system":"http://snomed.info/sct","code":"412739004","display":"BRCA2 gene mutation negative (finding)"},
              2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown (qualifier value)"}}
palb2_code = {"system":"http://loinc.org","code":"100761-6","display":"PALB2 gene targeted mutation analysis in Blood or Tissue by Molecular genetics method"}
chek2_code = {"system":"http://loinc.org","code":"89038-4","display":"CHEK2 gene deletion+duplication and full mutation analysis in Blood or Tissue by Molecular genetics method"}

#Method codes
pam50_code = {"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C120494","display":"PAM50"}

#Units and generic dictionaries
month_unit = {"value":"","unit":"month","system":"http://unitsofmeasure.org","code":"mo"}
unit_unit = {"value":"","unit":"Unit","system":"http://unitsofmeasure.org","code":"U"}
percent_unit = {"value":"","unit":"percent","system":"http://unitsofmeasure.org","code":"%"}
mgm2_unit = {"value":"","unit":"milligram per square meter","system":"http://unitsofmeasure.org","code":"mg/m2"}
mgkg_unit = {"value":"","unit":"milligram per kilogram","system":"http://unitsofmeasure.org","code":"mg/kg"}
mg_unit = {"value":"","unit":"milligram","system":"http://unitsofmeasure.org","code":"mg"}
auc2_unit = {"value":"","unit":"percent","system":"http://unitsofmeasure.org","code":""}
auc5_unit = {"value":"","unit":"percent","system":"http://unitsofmeasure.org","code":""}

chemo_unit_dict = {0:mgm2_unit,1:mgkg_unit,2:mg_unit,3:auc2_unit,4:auc5_unit,np.nan:unit_unit}

yesno_dict = {0:{"system":"http://snomed.info/sct","code":"373066001","display":"Yes"},
              1:{"system":"http://snomed.info/sct","code":"260415000","display":"Not detected"},
              2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
posnegun_dict = {0:{"system":"http://snomed.info/sct","code":"10828004","display":"Positive"},
                 1:{"system":"http://snomed.info/sct","code":"260385009","display":"Negative"},
                 2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
presence_dict = {0:{"system":"http://snomed.info/sct","code":"52101004","display":"Present (qualifier value)"},
                 1:{"system":"http://snomed.info/sct","code":"2667000","display":"Absent (qualifier value)"},
                 2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown (qualifier value)"}}
status_dict = {0: "not-done", 1: "stopped", 2: "completed", 3: "unknown"}
notrecorded = {"system":"http://snomed.info/sct", "code":"1220561009","display":"Not recorded (qualifier value)"}
#{"system":"","code":"","display":""}

# %%
#for this particular case, I need to create a function for lateralities
def lat(lat_index):
    laterality_dict = {0:{"system":"http://snomed.info/sct","code":"80248007","display":"Left breast structure (body structure)"},
                       1:{"system":"http://snomed.info/sct","code":"73056007","display":"Right breast structure (body structure)"}}
    notrecorded = {"system":"http://snomed.info/sct", "code":"1220561009","display":"Not recorded (qualifier value)"}
    if pd.isna(lat_index) == False and lat_index in laterality_dict:
        return laterality_dict[lat_index]
    else:
        return notrecorded
#define generic function for observation (valueCodeableConcept)
def obs(obs_index,obs_dict,obs_code,paturl,bundle,urnpat,urnint,urncase,category=np.nan,method=np.nan,bodysite=np.nan,effectivetimemonths=np.nan,effectivetimeweeks=np.nan,note=np.nan):
    obs_data = {}
    if pd.isna(obs_index) == False and obs_index in obs_dict:
        obs_data = {
            "status":"final",
            "code": {
                "coding":[obs_code]},
            "subject":{"reference":paturl},
            "valueCodeableConcept":{"coding":[obs_dict[obs_index]]}
            }
    #If not recorded, or giving error
    else:
        obs_data = {"status":"final",
                    "code":{
                        "coding":[obs_code]},
                    "subject":{"reference":paturl},
                    "valueCodeableConcept":{"coding":[notrecorded]}
                    }
        #Add category
    if pd.isna(category) == False:
        obs_data["category"] = {"coding":[category]}
        #Add method
    if pd.isna(method) == False:
        obs_data["method"] = {"coding":[method]}
        #Add bodysite
    if pd.isna(bodysite) == False:
        obs_data["bodySite"] = {"coding":[bodysite]}
        #Add month period
    if pd.isna(effectivetimemonths) == False and (isinstance(effectivetimemonths,int) == True or isinstance(effectivetimemonths,float) == True):
        obs_data["effectivePeriod"] = effectivemonths(effectivetimemonths)
        #Add week period
    if pd.isna(effectivetimeweeks) == False and (isinstance(effectivetimeweeks,int) == True or isinstance(effectivetimeweeks,float) == True):
        obs_data["effectivePeriod"] = effectiveweeks(effectivetimeweeks)
        #Add note
    if pd.isna(note) == False:
        obs_data["note"] = [{"text":note}]
    obs_resource = Observation(**obs_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(urnpat)+str(urnint)+str(urncase)).urn,"resource":obs_resource})

#define generic function for observation (valueQuantity)
def obsquant(obs_index,obs_unit,obs_code,paturl,bundle,urnpat,urnint,urncase,category=np.nan,method=np.nan,bodysite=np.nan,effectivetimemonths=np.nan,effectivetimeweeks=np.nan,note=np.nan,forceint=True):
    obs_data = {}
    obs_quantity_dict = obs_unit
    if pd.isna(obs_index) == False and obs_index != math.floor(obs_index) and (isinstance(obs_index,int) == True or isinstance(obs_index,float) == True):
        obs_quantity_dict["value"] = obs_index
        obs_data = {
            "status":"final",
            "code": {
                "coding":[obs_code]},
            "subject":{"reference":paturl},
            "valueQuantity":obs_quantity_dict
            }
    #if decimal is 0, transform into integer
    elif pd.isna(obs_index) == False and obs_index == math.floor(obs_index) and (isinstance(obs_index,int) == True or isinstance(obs_index,float) == True) and forceint == True:
        obs_quantity_dict["value"] = int(obs_index)
        obs_data = {
            "status":"final",
            "code": {
                "coding":[obs_code]},
            "subject":{"reference":paturl},
            "valueQuantity":obs_quantity_dict
            }
    #override the decimal to integer transformation, it stays decimal
    elif pd.isna(obs_index) == False and obs_index == math.floor(obs_index) and (isinstance(obs_index,int) == True or isinstance(obs_index,float) == True) and forceint == False:
        obs_quantity_dict["value"] = obs_index
        obs_data = {
            "status":"final",
            "code": {
                "coding":[obs_code]},
            "subject":{"reference":paturl},
            "valueQuantity":obs_quantity_dict
            }
    #If not recorded, or giving error
    else:
        obs_data = {"status":"final",
                    "code":{
                        "coding":[obs_code]},
                    "subject":{"reference":paturl},
                    "valueCodeableConcept":{"coding":[notrecorded]}
                    }
        #Add category
    if pd.isna(category) == False:
        obs_data["category"] = {"coding":[category]}
        #Add method
    if pd.isna(method) == False:
        obs_data["method"] = {"coding":[method]}
        #Add bodysite
    if pd.isna(bodysite) == False:
        obs_data["bodySite"] = {"coding":[bodysite]}
        #Add month period
    if pd.isna(effectivetimemonths) == False:
        obs_data["effectivePeriod"] = effectivemonths(effectivetimemonths)
        #Add week period
    if pd.isna(effectivetimeweeks) == False:
        obs_data["effectivePeriod"] = effectiveweeks(effectivetimeweeks)   
        #Add note
    if pd.isna(note) == False:
        obs_data["note"] = [{"text":note}]
    obs_resource = Observation(**obs_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(urnpat)+str(urnint)+str(urncase)).urn,"resource":obs_resource})

#define generic function for medication administration
def medadm(med_index,med_dict,paturl,bundle,urnpat,urnint,urncase,med_status=np.nan,med_status_dict=np.nan,med_dose=np.nan,med_unit=np.nan,med_unit_dict=np.nan,med_acc_dose=np.nan,med_acc_unit=np.nan,med_cycles=np.nan,effectivetimeweeks=np.nan,note=np.nan):
    medadm_data={}
    notrecorded = {"system":"http://snomed.info/sct", "code":"1220561009","display":"Not recorded (qualifier value)"}
    if pd.isna(med_status) == True or med_status not in med_status_dict:
        med_status = int(3)
    if pd.isna(med_status) == False and med_index in med_dict:
        medadm_data = {
            "status":med_status_dict[med_status],
            "medicationCodeableConcept":{"coding":[med_dict[med_index]]},
            "subject":{"reference":paturl}
        }
    else:
        medadm_data = {
            "status":med_status_dict[med_status],
            "medicationCodeableConcept":{"coding":[notrecorded]},
            "subject":{"reference":paturl}
        }
        #Add dose. Beware, if there are no units, the entire concept will be scrapped.
    if pd.isna(med_dose) == False and med_unit in med_unit_dict:
        dose_unit=med_unit
        dose_dict=med_unit_dict
        dose = dose_dict[dose_unit].copy()
        dose.update({"value":med_dose})
        medadm_data["dosage"] = {"dose":dose}
        #Add week period
    if pd.isna(effectivetimeweeks) == False and (isinstance(effectivetimeweeks,int) == True or isinstance(effectivetimeweeks,float) == True):
        medadm_data["effectivePeriod"] = effectiveweeks(effectivetimeweeks)
    else:
        medadm_data["effectivePeriod"] = effectiveweeks(0)
        #Prepare accummulated dose if available
    if pd.isna(med_acc_dose) == False and pd.isna(med_acc_unit) == False and med_acc_unit in med_unit_dict:
        acc_dose_unit = med_acc_unit
        acc_dose_dict = med_unit_dict
        acc_dose = acc_dose_dict[acc_dose_unit].copy()
        acc_dose["value"] = med_acc_dose
    else:
        acc_dose = np.nan
        #Add number of cycles and accummulated dose. If there is any issue with the dosage dictionary, it will override it.
            #Cycles == True and acc_dose and acc_units == True
    if (pd.isna(med_cycles) == False and (isinstance(med_cycles,int) or isinstance(med_cycles,float))) and (pd.isna(acc_dose) == False and pd.isna(med_acc_dose) == False and (isinstance(med_acc_dose,int) or isinstance(med_acc_dose,float))):
        try:
            medadm_data["dosage"]["rateRatio"] = {"numerator":acc_dose,"denominator":{"value":int(med_cycles),"unit":"Chemotherapy cycle","system":"http://snomed.info/sct","code":"399042005"}}
        except:
            medadm_data["dosage"] = {"rateRatio":{"numerator":acc_dose,"denominator":{"value":int(med_cycles),"unit":"Chemotherapy cycle","system":"http://snomed.info/sct","code":"399042005"}}}
            #Cycles == False and acc_dose and acc_units == True
    elif (pd.isna(med_cycles) == True or (isinstance(med_cycles,int) == False and isinstance(med_cycles,float)) == False) and (pd.isna(acc_dose) == False and pd.isna(med_acc_dose) == False and (isinstance(med_acc_dose,int) or isinstance(med_acc_dose,float))):
        try:
            medadm_data["dosage"]["rateRatio"] = {"numerator":acc_dose,"denominator":{"value":1,"unit":"Procedure","system":"http://snomed.info/sct","code":"71388002"}}
        except:
            medadm_data["dosage"] = {"rateRatio":{"numerator":acc_dose,"denominator":{"value":1,"unit":"Procedure","system":"http://snomed.info/sct","code":"71388002"}}}
            #Cycles == True and acc_dose and acc_units == False
    elif (pd.isna(med_cycles) == False and (isinstance(med_cycles,int) or isinstance(med_cycles,float))) and (pd.isna(acc_dose) == True or pd.isna(med_acc_dose) == True or (isinstance(med_acc_dose,int) == False and isinstance(med_acc_dose,float) == False)):
        try:
            medadm_data["dosage"]["rateQuantity"] = {"value":int(med_cycles),"unit":"Chemotherapy cycle","system":"http://snomed.info/sct","code":"399042005"}
        except:
            medadm_data["dosage"] = {"rateQuantity":{"value":int(med_cycles),"unit":"Chemotherapy cycle","system":"http://snomed.info/sct","code":"399042005"}}   
        #Add note
    if pd.isna(note) == False:
        medadm_data["note"] = [{"text":note}]
    medadm_resource = MedicationAdministration(**medadm_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(urnpat)+str(urnint)+str(urncase)).urn,"resource":medadm_resource})

#define function for histopat
def hist(hist_index,paturl,bundle,urnpat,urnint,urncase,bodysite=np.nan,note=np.nan):
    hist_data = {}
    if pd.isna(hist_index) == False and hist_index in histopat_dict:
        hist_data = {
            "code":{"coding":[histopat_dict[hist_index]]},
            "subject":{"reference":paturl}
        }
    else:
        hist_data = {
            "code":{"coding":[notrecorded]},
            "subject":{"reference":paturl}
        }
    #Add bodysite. For some reason, in condition, bodysite has to be a list containing the coding dictionary, instead of the dictionary as a whole
    if pd.isna(bodysite) == False:
        hist_data["bodySite"] = [{"coding":[bodysite]}]
    #add note
    if pd.isna(note) == False:
        hist_data["note"] = [{"text":note}]

    histopat = Condition(**hist_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(urnpat)+str(urnint)+str(urncase)).urn,"resource":histopat})
        
#define function for effectivetimePeriod in weeks
def effectiveweeks(weeks):
    startdate = np.datetime64('1970-01-01')
    startdatetext = str(startdate)+'T00:00:00+00:00'
    try:
        timepassed = np.timedelta64(weeks,'W')
    except:
        if pd.isna(weeks) == False and isinstance(weeks,float) == True and weeks == math.floor(weeks):
           timepassed = np.timedelta64(int(weeks),'W') 
        elif pd.isna(weeks) == False and isinstance(weeks,float) == True and weeks != math.floor(weeks):
            days = math.ceil(weeks*7)
            timepassed = np.timedelta64(days,'D')
        else:
            timepassed = np.timedelta64(0,'W')
    finaldate = startdate+timepassed
    finaldatetext = str(finaldate)+'T00:00:00+00:00'
    effectivePeriodweeks = {"start":startdatetext,
                            "end":finaldatetext}
    return effectivePeriodweeks

#define function for effectivetimePeriod in months
def effectivemonths(months):
    startdate = np.datetime64('1970-01')
    startdatetext = str(startdate)+'-01T00:00:00+00:00'
    try:
        timepassed = np.timedelta64(months,'M')
        finaldate = startdate+timepassed
        finaldatetext = str(finaldate)+'-01T00:00:00+00:00'
    except:
        if pd.isna(months) == False and isinstance(months,float) == True and months == math.floor(months):
            timepassed = np.timedelta64(int(months),'M')
            finaldate = startdate+timepassed
            finaldatetext = str(finaldate)+'-01T00:00:00+00:00'
        elif pd.isna(months) == False and isinstance(months,float) == True and months != math.floor(months):
            weeks = math.ceil(months*4.34524)
            startdate = np.datetime64('1970-01-01')
            timepassed = np.timedelta64(weeks,'W')
            finaldate = startdate+timepassed
            finaldatetext = str(finaldate)+'T00:00:00+00:00'
        else:
            timepassed = np.timedelta64(0,'M')
            finaldate = startdate+timepassed
            finaldatetext = str(finaldate)+'-01T00:00:00+00:00'
    effectivePeriodmonths = {"start":startdatetext,
                             "end":finaldatetext}
    return effectivePeriodmonths

# %%
for index in range(len(data_csv)):
    #Empty variables not to reuse older data
    d_index = {}
    #for i in range(len(d_csv)):
    #    d_index[data_csv.columns[i]+"_index"] = np.nan
    
    #Data taken from lists
    for i in range(len(d_csv)):
        d_index[data_csv.columns[i]+"_index"] = d_csv[data_csv.columns[i]+"_csv"][index]
        
    #Bundle setup
    bundle_data = {
        "id": str("UC68-bundle-")+str(d_index["redcap_data_access_group_index"])+"-"+str(index+1),
        "type": "collection",
        "entry":[],
    }
    bundle = Bundle(**bundle_data)
    
    #Patient setup
    paturl = uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("patient")+str(68)).urn
    pat_data = {
        "identifier":[
            {
            "value": d_index["record_id_index"],
            "use": "secondary"
            }
        ],
    }
        
    pat= Patient(**pat_data)
    bundle.entry.append({"fullUrl":paturl,"resource":pat})

    #sex setup
    obs(d_index["sex_index"],sex_dict,sex_code,paturl,bundle,d_index["record_id_index"],"sex",68)
    
    #primary condition setup
    primary_condition_data ={
        "code":{"coding":[patientclass_dict[d_index["patientclass_index"]]]},
        "subject":{"reference":paturl},
        "bodySite":[
            {
                "coding":[
                    {
                        "system":"http://snomed.info/sct",
                        "code":"76752008",
                        "display":"Breast structure (body structure)"
                    }
                ]
            }
        ]
    }
    #onsetage setup
    if pd.isna(d_index["onsetage_index"]) == False and (isinstance(d_index["onsetage_index"],int) == True or isinstance(d_index["onsetage_index"],float) == True):
        primary_condition_data["onsetAge"]= {
                    "value": d_index["onsetage_index"],
                    "unit": "years",
                    "code": "a",
                    "system": "http://unitsofmeasure.org"}
    primary_condition = Condition(**primary_condition_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("condition")+str(68)).urn,"resource":primary_condition})
    #menop setup
    obs(d_index["menop_index"],menop_dict,menop_code,paturl,bundle,d_index["record_id_index"],"menop",68)
    #n_preg setup
    obsquant(d_index["n_preg_index"],unit_unit,n_preg_code,paturl,bundle,d_index["record_id_index"],"n_preg",68)
    #lactation and lactation_dur setup
    obs(d_index["lactation_index"],yesno_dict,lactation_code,paturl,bundle,d_index["record_id_index"],"lactation",68,effectivetimemonths=d_index["lactation_dur_index"],note="If used in this observation, the effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates.")
    #symptoms setup
    obs(d_index["symptoms_index"],yesno_dict,symptoms_code,paturl,bundle,d_index["record_id_index"],"symptoms",68)
    #screening setup
    obs(d_index["screening_index"],yesno_dict,screening_code,paturl,bundle,d_index["record_id_index"],"screening",68)
    #famhisto_b setup
    obs(d_index["famhisto_b_index"],famhisto_dict,famhisto_b_code,paturl,bundle,d_index["record_id_index"],"famhisto_b",68)
    #famhisto_o setup
    obs(d_index["famhisto_o_index"],famhisto_dict,famhisto_o_code,paturl,bundle,d_index["record_id_index"],"famhisto_o",68)
    #Laterality 1
    #bcstproxy setup
    bcstproxy_data = {}
        #triple negative
    if pd.isna(d_index["bcstproxy_index"]) == False and d_index["bcstproxy_index"] == 0:
        bcstproxy_data = {
            "status":"final",
            "code":{"coding":[bcstproxy_code]},
            "subject":{"reference":paturl},
            "valueCodeableConcept":{"coding":[bc_dict[0]]}
        }
        #non triple negative
    elif d_index["bcstproxy_index"] in bc_dict and d_index["bcstproxy_index"] > 0:
        bcstproxy_data = {
            "status":"final",
            "code":{"coding":[bcstproxy_code]},
            "subject":{"reference":paturl},
            "component":bc_dict[d_index["bcstproxy_index"]]
        }
        #not recorded
    else:
        bcstproxy_data = {
            "status":"final",
            "code":{"coding":[bcstproxy_code]},
            "subject":{"reference":paturl},
            "valueCodeableConcept":{"coding":[notrecorded]}
        }
    if pd.isna(d_index["bcstproxy_index"]) == False and d_index["bcstproxy_index"] in bc_note_dict:
        bcstproxy_data["note"] = [{"text":bc_note_dict[d_index["bcstproxy_index"]]}]
        #laterality
    bcstproxy_data["bodySite"] = {"coding":[lat(d_index["laterality_index"])]}
    bctsproxy = Observation(**bcstproxy_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("bcstproxy")+str(68)).urn,"resource":bctsproxy})
    #bcst_pam50 setup
    bcst_pam50_data = {}
        #triple negative
    if pd.isna(d_index["bcst_pam50_index"]) == False and d_index["bcst_pam50_index"] == 0:
        bcst_pam50_data = {
            "status":"final",
            "code":{"coding":[bcst_pam50_code]},
            "subject":{"reference":paturl},
            "method": {"coding":[pam50_code]},
            "valueCodeableConcept":{"coding":[bc_dict[0]]}
        }
        #non triple negative
    elif d_index["bcst_pam50_index"] in bc_dict and d_index["bcst_pam50_index"] > 0:
        bcst_pam50_data = {
            "status":"final",
            "code":{"coding":[bcst_pam50_code]},
            "subject":{"reference":paturl},
            "method": {"coding":[pam50_code]},
            "component":bc_dict[d_index["bcst_pam50_index"]]
        }
        #not recorded
    else:
        bcst_pam50_data = {
            "status":"final",
            "code":{"coding":[bcst_pam50_code]},
            "subject":{"reference":paturl},
            "method": {"coding":[pam50_code]},
            "valueCodeableConcept":{"coding":[notrecorded]}
        }
    if pd.isna(d_index["bcst_pam50_index"]) == False and d_index["bcst_pam50_index"] in bc_note_dict:
        bcst_pam50_data["note"] = [{"text":bc_note_dict[d_index["bcst_pam50_index"]]}]
        #laterality
    bcst_pam50_data["bodySite"] = {"coding":[lat(d_index["laterality_index"])]}
    bcst_pam50 = Observation(**bcst_pam50_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("bcst_pam50")+str(68)).urn,"resource":bcst_pam50})
    #histopat 1 setup
    hist(d_index["histopat_index"],paturl,bundle,d_index["record_id_index"],"histopat1",68,bodysite=lat(d_index["laterality_index"]),note="Histopathological type")
    #histopat 2 setup
    hist(d_index["histopat_2_index"],paturl,bundle,d_index["record_id_index"],"histopat2",68,bodysite=lat(d_index["laterality_index"]),note="Additional histopathological type")
    #histopat 3 setup
    hist(d_index["histopat_3_index"],paturl,bundle,d_index["record_id_index"],"histopat3",68,bodysite=lat(d_index["laterality_index"]),note="Second additional Histopathological type")
    #tnm_ct setup
    obs(d_index["tnm_ct_index"],tnm_ct_dict,tnm_ct_code,paturl,bundle,d_index["record_id_index"],"tnm_ct",68,bodysite=lat(d_index["laterality_index"]))
    #tnm_cn setup
    obs(d_index["tnm_cn_index"],tnm_cn_dict,tnm_cn_code,paturl,bundle,d_index["record_id_index"],"tnm_cn",68,bodysite=lat(d_index["laterality_index"]))
    #tnm_cm setup
    obs(d_index["tnm_cm_index"],tnm_cm_dict,tnm_cm_code,paturl,bundle,d_index["record_id_index"],"tnm_cm",68,bodysite=lat(d_index["laterality_index"]))
    #grade setup
    obs(d_index["grade_index"],grade_dict,grade_code,paturl,bundle,d_index["record_id_index"],"grade",68,bodysite=lat(d_index["laterality_index"]))
    #dcis setup
    obs(d_index["dcis_index"],dcis_dict,dcis_code,paturl,bundle,d_index["record_id_index"],"dcis",68,bodysite=lat(d_index["laterality_index"]))
    #er setup
    obsquant(d_index["er_index"],percent_unit,er_code,paturl,bundle,d_index["record_id_index"],"er",68,bodysite=lat(d_index["laterality_index"]),forceint=False)
    #pr setup
    obsquant(d_index["pr_index"],percent_unit,pr_code,paturl,bundle,d_index["record_id_index"],"pr",68,bodysite=lat(d_index["laterality_index"]),forceint=False)
    #her2ihc setup
    obs(d_index["her2ihc_index"],her2ihc_dict,her2ihc_code,paturl,bundle,d_index["record_id_index"],"her2ihc",68,bodysite=lat(d_index["laterality_index"]))
    #ki67 setup
    obsquant(d_index["ki67_index"],percent_unit,ki67_code,paturl,bundle,d_index["record_id_index"],"ki67",68,bodysite=lat(d_index["laterality_index"]),forceint=False)
    #her2fish setup
    obs(d_index["her2fish_index"],posnegun_dict,her2fish_code,paturl,bundle,d_index["record_id_index"],"her2fish",68,bodysite=lat(d_index["laterality_index"]))
    if pd.isna(d_index["laterality_2_index"]) == False and d_index["laterality_2_index"] in laterality_dict:
    #if True == True:
        #bcstproxy_2 setup
        bcstproxy_2_data = {}
            #triple negative
        if pd.isna(d_index["bcstproxy_2_index"]) == False and d_index["bcstproxy_2_index"] == 0:
            bcstproxy_2_data = {
                "status":"final",
                "code":{"coding":[bcstproxy_code]},
                "subject":{"reference":paturl},
                "valueCodeableConcept":{"coding":[bc_dict[0]]}
            }
            #non triple negative
        elif d_index["bcstproxy_2_index"] in bc_dict and d_index["bcstproxy_2_index"] > 0:
            bcstproxy_2_data = {
                "status":"final",
                "code":{"coding":[bcstproxy_code]},
                "subject":{"reference":paturl},
                "component":bc_dict[d_index["bcstproxy_2_index"]]
            }
            #not recorded
        else:
            bcstproxy_2_data = {
                "status":"final",
                "code":{"coding":[bcstproxy_code]},
                "subject":{"reference":paturl},
                "valueCodeableConcept":{"coding":[notrecorded]}
            }
        if pd.isna(d_index["bcstproxy_2_index"]) == False and d_index["bcstproxy_2_index"] in bc_note_dict:
            bcstproxy_2_data["note"] = [{"text":bc_note_dict[d_index["bcstproxy_2_index"]]}]
            #laterality
        bcstproxy_2_data["bodySite"] = {"coding":[lat(d_index["laterality_2_index"])]}
        bctsproxy_2 = Observation(**bcstproxy_2_data)
        bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("bcstproxy_2")+str(68)).urn,"resource":bctsproxy_2})
        #bcst_pam50_2 setup
        bcst_pam50_2_data = {}
            #triple negative
        if pd.isna(d_index["bcst_pam50_2_index"]) == False and d_index["bcst_pam50_2_index"] == 0:
            bcst_pam50_2_data = {
                "status":"final",
                "code":{"coding":[bcst_pam50_code]},
                "subject":{"reference":paturl},
                "method": {"coding":[pam50_code]},
                "valueCodeableConcept":{"coding":[bc_dict[0]]}
            }
            #non triple negative
        elif d_index["bcst_pam50_2_index"] in bc_dict and d_index["bcst_pam50_2_index"] > 0:
            bcst_pam50_2_data = {
                "status":"final",
                "code":{"coding":[bcst_pam50_code]},
                "subject":{"reference":paturl},
                "method": {"coding":[pam50_code]},
                "component":bc_dict[d_index["bcst_pam50_2_index"]]
            }
            #not recorded
        else:
            bcst_pam50_2_data = {
                "status":"final",
                "code":{"coding":[bcst_pam50_code]},
                "subject":{"reference":paturl},
                "method": {"coding":[pam50_code]},
                "valueCodeableConcept":{"coding":[notrecorded]}
            }
        if pd.isna(d_index["bcst_pam50_2_index"]) == False and d_index["bcst_pam50_2_index"] in bc_note_dict:
            bcst_pam50_2_data["note"] = [{"text":bc_note_dict[d_index["bcst_pam50_2_index"]]}]
            #laterality
        bcst_pam50_2_data["bodySite"] = {"coding":[lat(d_index["laterality_2_index"])]}
        bcst_pam50_2 = Observation(**bcst_pam50_2_data)
        bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("bcst_pam50_2")+str(68)).urn,"resource":bcst_pam50_2})
        #histopat_1_2 setup
        hist(d_index["histopat_1_2_index"],paturl,bundle,d_index["record_id_index"],"histopat1_2",68,bodysite=lat(d_index["laterality_2_index"]),note="Histopathological type")
        #histopat_2_2 setup
        hist(d_index["histopat_2_2_index"],paturl,bundle,d_index["record_id_index"],"histopat2_2",68,bodysite=lat(d_index["laterality_2_index"]),note="Additional histopathological type")
        #histopat_3_2 setup
        hist(d_index["histopat_3_2_index"],paturl,bundle,d_index["record_id_index"],"histopat3_2",68,bodysite=lat(d_index["laterality_2_index"]),note="Second additional Histopathological type")
        #tnm_ct_2 setup
        obs(d_index["tnm_ct_2_index"],tnm_ct_dict,tnm_ct_code,paturl,bundle,d_index["record_id_index"],"tnm_ct_2",68,bodysite=lat(d_index["laterality_2_index"]))
        #tnm_cn_2 setup
        obs(d_index["tnm_cn_2_index"],tnm_cn_dict,tnm_cn_code,paturl,bundle,d_index["record_id_index"],"tnm_cn_2",68,bodysite=lat(d_index["laterality_2_index"]))
        #tnm_cm_2 setup
        obs(d_index["tnm_cm_2_index"],tnm_cm_dict,tnm_cm_code,paturl,bundle,d_index["record_id_index"],"tnm_cm_2",68,bodysite=lat(d_index["laterality_2_index"]))
        #grade_2 setup
        obs(d_index["grade_2_index"],grade_dict,grade_code,paturl,bundle,d_index["record_id_index"],"grade_2",68,bodysite=lat(d_index["laterality_2_index"]))
        #dcis_2 setup
        obs(d_index["dcis_2_index"],dcis_dict,dcis_code,paturl,bundle,d_index["record_id_index"],"dcis_2",68,bodysite=lat(d_index["laterality_2_index"]))
        #er_2 setup
        obsquant(d_index["er_2_index"],percent_unit,er_code,paturl,bundle,d_index["record_id_index"],"er_2",68,bodysite=lat(d_index["laterality_2_index"]),forceint=False)
        #pr_2 setup
        obsquant(d_index["pr_2_index"],percent_unit,pr_code,paturl,bundle,d_index["record_id_index"],"pr_2",68,bodysite=lat(d_index["laterality_2_index"]),forceint=False)
        #her2ihc_2 setup
        obs(d_index["her2ihc_2_index"],her2ihc_dict,her2ihc_code,paturl,bundle,d_index["record_id_index"],"her2ihc_2",68,bodysite=lat(d_index["laterality_2_index"]))
        #ki67_2 setup
        obsquant(d_index["ki67_2_index"],percent_unit,ki67_code,paturl,bundle,d_index["record_id_index"],"ki67_2",68,bodysite=lat(d_index["laterality_2_index"]),forceint=False)
        #her2fish_2 setup
        obs(d_index["her2fish_2_index"],posnegun_dict,her2fish_code,paturl,bundle,d_index["record_id_index"],"her2fish_2",68,bodysite=lat(d_index["laterality_2_index"]))
    #hcontrac and hcohtrac_dur setup
    obs(d_index["hcontrac_index"],yesno_dict,hcontrac_code,paturl,bundle,d_index["record_id_index"],"hcontrac",68,effectivetimemonths=d_index["hcohtrac_dur_index"],note="If used in this observation, the effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates.")
    #hormtherapy and hormtherapy_dur setup
    obs(d_index["hormtherapy_index"],yesno_dict,hormtherapy_code,paturl,bundle,d_index["record_id_index"],"hormtherapy",68,effectivetimemonths=d_index["hormtherapy_dur_index"],note="If used in this observation, the effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates.")
    #brca1 setup
    obs(d_index["brca1_index"],posnegun_dict,brca1_code,paturl,bundle,d_index["record_id_index"],"brca1",68)
    #brca2 setup
    obs(d_index["brca2_index"],posnegun_dict,brca2_code,paturl,bundle,d_index["record_id_index"],"brca2",68)
    #palb2 setup
    obs(d_index["palb2_index"],posnegun_dict,palb2_code,paturl,bundle,d_index["record_id_index"],"palb2",68)
    #chek2 setup
    obs(d_index["chek2_index"],posnegun_dict,chek2_code,paturl,bundle,d_index["record_id_index"],"chek2",68)
    
    #Making JSON readable and parsing it
    prettybundle = json.loads(bundle.json())
    parsed = json.dumps(prettybundle, indent=4)
    #Giving JSON an export name and path
    relative_path_json = "/UC6&8/OUTPUT/"
    if os.path.exists(absolute_path+relative_path_json) == False:
        os.makedirs(absolute_path+relative_path_json)
    if pd.isna(d_index["redcap_data_access_group_index"]) == False:
        dag_path_json = str(d_index["redcap_data_access_group_index"])+"/"
        if os.path.exists(absolute_path+relative_path_json+dag_path_json) == False:
            os.makedirs(absolute_path+relative_path_json+dag_path_json)
    else:
        dag_path_json = ""
    jsonfilename = str(d_index["record_id_index"])+".json"
    jsonpath = absolute_path+relative_path_json+dag_path_json+jsonfilename
    #Writing JSON
    with open(jsonpath, "w") as outfile:
        outfile.write(parsed)

