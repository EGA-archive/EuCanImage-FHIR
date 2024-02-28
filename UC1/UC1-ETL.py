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
import json

# %%
#Import data from CSV. Change formats to readable ones by dictionaries.
#absolute_path = os.path.dirname(os.getcwd())
absolute_path = os.getcwd()
relative_path_csv = "/UC1/CSV/UseCase1_testdata.csv"
csvpath = absolute_path+relative_path_csv
data_csv = pd.read_csv(csvpath)
data_csv = data_csv.astype({
    "sex":"Int64",
    "indeterm_nodule":"Int64",
    "diag_method":"Int64",
    "hep_b":"Int64",
    "hep_c":"Int64",
    "nafld":"Int64",
    "ald":"Int64",
    "autoimm_hep":"Int64",
    "psc":"Int64",
    "others":"Int64",
    "chronic_hep":"Int64",
    "cirrhosis":"Int64",
    "child_pugh": "Int64"
})
#Convert table to dictionary containing all variables as lists
d_csv = {}
for i in range(len(data_csv.columns)):
    d_csv[data_csv.columns[i]+"_csv"] = data_csv[data_csv.columns[i]].tolist()

# %%
#Data Dictionaries
sex_code= {"system": "http://loinc.org","code": "76689-9","display": "Sex assigned at birth"}
sex_dict = {0: {"system":"http://loinc.org","code":"LA2-8","display":"Male"},
            1: {"system":"http://loinc.org","code":"LA3-6","display":"Female"}}
indeterm_nodule_dict = {0:{"system":"http://terminology.hl7.org/CodeSystem/condition-ver-status","code":"confirmed","display":"Confirmed"},
                        1:{"system":"http://terminology.hl7.org/CodeSystem/condition-ver-status","code":"refuted","display":"Refuted"}}
diag_method_dict = {0:{"system":"http://snomed.info/sct","code":"394597005","display":"Histopathology (qualifier value)"},
                    1:{"system":"http://snomed.info/sct","code":"360037004","display":"Imaging - action (qualifier value)"}}
hep_b_code = {"system":"http://snomed.info/sct","code":"365863005","display":"Finding of Hepatitis B status (finding)"}
hep_c_code = {"system":"http://snomed.info/sct","code":"365865003","display":"Finding of Hepatitis C status (finding)"}
ald_code = {"system":"http://snomed.info/sct","code":"41309000","display":"Alcoholic liver damage (disorder)"}
autoimm_hep_code = {"system":"http://snomed.info/sct","code":"408335007","display":"Autoimmune hepatitis (disorder)"}
nafld_code = {"system":"http://snomed.info/sct","code":"1231824009","display":"Non-alcoholic fatty liver disease (disorder)"}
psc_code = {"system":"http://snomed.info/sct","code":"197441003","display":"Primary sclerosing cholangitis (disorder)"}
others_code = {"system":"http://snomed.info/sct","code":"235856003:42752001=74964007","display":"Hepatopathy where Due to = Other"}
chronic_hep_code = {"system":"http://snomed.info/sct","code":"76783007","display":"Chronic hepatitis (disorder)"}
cirrhosis_code = {"system":"http://snomed.info/sct","code":"19943007","display":"Cirrhosis of liver (disorder)"}
child_pugh_code = {"system":"http://snomed.info/sct","code":"3191000175106","display":"Child-Pugh score (assessment scale)"}
child_pugh_dict = {0:{"system":"http://snomed.info/sct","code":"710065009","display":"Child-Pugh score class A (finding)"},
                   1:{"system":"http://snomed.info/sct","code":"710066005","display":"Child-Pugh score class B (finding)"},
                   2:{"system":"http://snomed.info/sct","code":"710067001","display":"Child-Pugh score class C (finding)"}}
posnegun_dict = {0:{"system":"http://snomed.info/sct","code":"10828004","display":"Positive"},
       1:{"system":"http://snomed.info/sct","code":"260385009","display":"Negative"},
       2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
yesno_dict = {0:{"system":"http://snomed.info/sct","code":"373066001","display":"Yes"},
              1:{"system":"http://snomed.info/sct","code":"260415000","display":"Not detected"},
              2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
notrecorded = {"system":"http://snomed.info/sct", "code":"1220561009","display":"Not recorded (qualifier value)"}
#{"system":"","code":"","display":""}

# %%
#define generic function for observation (valueCodeableConcept)
def obs(obs_index,obs_dict,obs_code,paturl,bundle,urnpat,urnint,urncase,category=np.nan,method=np.nan,evidence=np.nan,bodysite=np.nan,effectivetimemonths=np.nan,effectivetimeweeks=np.nan,note=np.nan):
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
        #Add evidence
    if pd.isna(evidence) == False:
        obs_data["evidence"] = {"coding":[evidence]}
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
def obsquant(obs_index,obs_unit,obs_code,paturl,bundle,urnpat,urnint,urncase,category=np.nan,method=np.nan,evidence=np.nan,bodysite=np.nan,effectivetimemonths=np.nan,effectivetimeweeks=np.nan,note=np.nan,forceint=True):
    obs_data = {}
    obs_quantity_dict = obs_unit.copy()
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
        #Add evidence
    if pd.isna(evidence) == False:
        obs_data["evidence"]= {"coding":[evidence]}
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
        "id": str("UC1-bundle-")+str(d_index["redcap_data_access_group_index"])+"-"+str(index+1),
        "type": "collection",
        "entry":[],
    }
    bundle = Bundle(**bundle_data)
    
    #Patient setup
    paturl = uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("patient")+str(1)).urn
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
    obs(d_index["sex_index"],sex_dict,sex_code,paturl,bundle,d_index["record_id_index"],"sex",1)
    
    #primary condition setup
    primary_condition_data ={
        "code":{"coding":[{
                    "system":"http://snomed.info/sct",
                    "code":"95214007",
                    "display":"Primary malignant neoplasm of liver (disorder)"
        }]},
        "subject":{"reference":paturl},
        "bodySite":[
            {
                "coding":[
                    {
                        "system":"http://snomed.info/sct",
                        "code":"10200004",
                        "display":"Liver structure (body structure)"
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
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("condition")+str(1)).urn,"resource":primary_condition})
    
    #indeterm_nodule setup
    indeterm_nodule_data ={
        "code":{"coding":[{
                    "system":"http://snomed.info/sct",
                    "code":"109841003",
                    "display":"Liver cell carcinoma (disorder)"
        }]},
        "subject":{"reference":paturl}}
        #variable setup
    if pd.isna(d_index["indeterm_nodule_index"]) == False and d_index["indeterm_nodule_index"] in indeterm_nodule_dict:
        indeterm_nodule_data["verificationStatus"] = {"coding":[indeterm_nodule_dict[d_index["indeterm_nodule_index"]]]}
    else:
        indeterm_nodule_data["verificationStatus"] = {"coding":[{"system":"http://terminology.hl7.org/CodeSystem/condition-ver-status","code":"unconfirmed","display":"Unconfirmed"}]}
        #evidence setup
    if pd.isna(d_index["diag_method_index"]) == False and d_index["diag_method_index"] in diag_method_dict:
        indeterm_nodule_data["evidence"] = [{"code":[{"coding":[diag_method_dict[d_index["diag_method_index"]]]}]}]
    else:
        indeterm_nodule_data["evidence"] = [{"code":[{"coding":[notrecorded]}]}]
    
    indeterm_nodule = Condition(**indeterm_nodule_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("indeterm_nodule")+str(1)).urn,"resource":indeterm_nodule})

    #Hep B setup
    obs(d_index["hep_b_index"],posnegun_dict,hep_b_code,paturl,bundle,d_index["record_id_index"],"hep_b",1)
    #Hep C setup
    obs(d_index["hep_c_index"],posnegun_dict,hep_c_code,paturl,bundle,d_index["record_id_index"],"hep_c",1)
    #ald setup
    obs(d_index["ald_index"],yesno_dict,ald_code,paturl,bundle,d_index["record_id_index"],"ald",1)
    #autoimm_hep setup
    obs(d_index["autoimm_hep_index"],yesno_dict,autoimm_hep_code,paturl,bundle,d_index["record_id_index"],"autoimm_hep",1)
    #nafld setup
    obs(d_index["nafld_index"],yesno_dict,nafld_code,paturl,bundle,d_index["record_id_index"],"nafld",1)
    #psc setup
    obs(d_index["psc_index"],yesno_dict,nafld_code,paturl,bundle,d_index["record_id_index"],"psc",1)
    #others setup
    obs(d_index["others_index"],yesno_dict,others_code,paturl,bundle,d_index["record_id_index"],"others",1)
    #chronic_hep setup
    obs(d_index["chronic_hep_index"],yesno_dict,chronic_hep_code,paturl,bundle,d_index["record_id_index"],"chronic_hep",1)
    #cirrhosis setup
    obs(d_index["cirrhosis_index"],yesno_dict,cirrhosis_code,paturl,bundle,d_index["record_id_index"],"cirrhosis",1)
    #child_pugh setup
    obs(d_index["child_pugh_index"],child_pugh_dict,child_pugh_code,paturl,bundle,d_index["record_id_index"],"child_pugh",1)
            
    #Making JSON readable and parsing it
    prettybundle = json.loads(bundle.json())
    parsed = json.dumps(prettybundle, indent=4)
    #Giving JSON an export name and path
    relative_path_json = "/UC1/OUTPUT/"
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