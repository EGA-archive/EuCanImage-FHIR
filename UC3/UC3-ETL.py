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
relative_path_csv = "/UC3/CSV/UseCase3_testdata.csv"
csvpath = absolute_path+relative_path_csv
data_csv = pd.read_csv(csvpath)
data_csv = data_csv.astype({
    "sex":"Int64",
    "histopat":"Int64",
    "diag_method":"Int64",
    "comorb_present":"Int64",
    "bca":"Int64",
    "biloma":"Int64",
    "cca":"Int64",
    "dn":"Int64",
    "fhcc":"Int64",
    "fhs":"Int64",
    "fnh":"Int64",
    "hemangioma":"Int64",
    "abscess":"Int64",
    "adenoma":"Int64",
    "pha":"Int64",
    "hhcyst":"Int64",
    "hcc":"Int64",
    "cirrhosis":"Int64",
    "phl":"Int64",
    "shcyst":"Int64",
    "thad":"Int64",
    "tnm_pt":"Int64",
    "tnm_ct":"Int64",
    "tnm_pn":"Int64",
    "tnm_cn":"Int64",
    "tnm_pm":"Int64",
    "tnm_cm":"Int64"
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
histopat_dict = {0:{"system":"http://snomed.info/sct","code":"94381002:260716003=93761005","display":"Secondary malignant neoplasm of liver where Site of primary tumour = Primary malignant neoplasm of colon"},
                 1:{"system":"http://snomed.info/sct","code":"52988006:410657003=74964007","display":"Lesion where Type - attribute = Other"},
                 2:{"system":"http://snomed.info/sct","code":"52988006:410657003=260415000","display":"Lesion where Type - = Not detected"}}
diag_method_dict = {0:{"system":"http://snomed.info/sct","code":"394597005","display":"Histopathology (qualifier value)"},
                    1:{"system":"http://snomed.info/sct","code":"360037004","display":"Imaging - action (qualifier value)"}}
comorb_code = {"system":"http://snomed.info/sct","code":"398192003","display":"Co-morbid conditions (finding)"}
bca_code = {"system":"http://snomed.info/sct", "code":"448997008","display":"Cystadenoma of liver (disorder)"}
biloma_code = {"system":"http://snomed.info/sct","code":"445137006","display":"Intraabdominal bile collection"}
cca_code = {"system":"http://snomed.info/sct","code":"70179006","display":"Cholangiocarcinoma (morphologic abnormality)"}
dn_code = {"system":"http://snomed.info/sct","code":"448149008","display":"Dysplastic nodule (morphologic abnormality)"}
fhcc_code = {"system":"http://snomed.info/sct","code":"253018005","display":"Fibrolamellar hepatocellular carcinoma (disorder)"}
fhs_code = {"system":"http://snomed.info/sct","code":"197321007:116676008=860692001","display":"Steatosis of liver where Associated morphology = Focal damage"}
fnh_code = {"system":"http://snomed.info/sct","code":"278527001","display":"Focal nodular hyperplasia of liver"}
hemangioma_code = {"system":"http://snomed.info/sct","code":"400210000","display":"Hemangioma (disorder)"}
abscess_code = {"system":"http://snomed.info/sct","code":"27916005","display":"Abscess of liver (disorder)"}
adenoma_code = {"system":"http://snomed.info/sct","code":"424263008","display":"Adenoma of liver (disorder)"}
pha_code = {"system":"http://snomed.info/sct","code":"109844006","display":"Angiosarcoma of liver (disorder)"}
hhcyst_code = {"system":"http://snomed.info/sct","code":"26103000","display":"Echinococcosis of liver (disorder)"}
cirrhosis_code = {"system":"http://snomed.info/sct","code":"19943007","display":"Cirrhosis of liver (disorder)"}
phl_code = {"system":"http://snomed.info/sct","code":"1153383006","display":"Malignant lymphoma of liver (disorder)"}
shcyst_code = {"system":"http://snomed.info/sct","code":"85057007","display":"Cyst of liver (disorder)"}
#PLACEHOLDER FOR THAD
thad_code = {"system":"http://snomed.info/sct","code":"398192003","display":"Co-morbid conditions"}
hcc_code = {"system":"http://snomed.info/sct","code":"109841003","display":"Liver cell carcinoma (disorder)"}
tnm_pt_code = {"system":"http://snomed.info/sct","code":"1222589003","display":"American Joint Committee on Cancer pathological T category allowable value (qualifier value)"}
tnm_pt_dict = {0:{"system":"http://snomed.info/sct","code":"1228950008","display":"pTX"},
               1:{"system":"http://snomed.info/sct","code":"1228951007","display":"pT0"},
               2:{"system":"http://snomed.info/sct","code":"1228953005","display":"pTis"},
               3:{"system":"http://snomed.info/sct","code":"1228957006","display":"pT1"},
               4:{"system":"http://snomed.info/sct","code":"1229852009","display":"pT2"},
               5:{"system":"http://snomed.info/sct","code":"1229859000","display":"pT3"},
               6:{"system":"http://snomed.info/sct","code":"1229864001","display":"pT4"},
               7:{"system":"http://snomed.info/sct","code":"1229865000","display":"pT4a"},
               8:{"system":"http://snomed.info/sct","code":"1229866004","display":"pT4b"}}
tnm_ct_code = {"system":"http://snomed.info/sct","code":"1222585009","display":"American Joint Committee on Cancer clinical T category allowable value (qualifier value)"}
tnm_ct_dict = {0:{"system":"http://snomed.info/sct","code":"1222604002","display":"cTX"},
               1:{"system":"http://snomed.info/sct","code":"1228882005","display":"cT0"},
               2:{"system":"http://snomed.info/sct","code":"1228884006","display":"cTis"},
               3:{"system":"http://snomed.info/sct","code":"1228889001","display":"cT1"},
               4:{"system":"http://snomed.info/sct","code":"1228929004","display":"cT2"},
               5:{"system":"http://snomed.info/sct","code":"1228938002","display":"cT3"},
               6:{"system":"http://snomed.info/sct","code":"1228944003","display":"cT4"},
               7:{"system":"http://snomed.info/sct","code":"1228945002","display":"cT4a"},
               8:{"system":"http://snomed.info/sct","code":"1228946001","display":"cT4b"}}
tnm_pn_code = {"system":"http://snomed.info/sct","code":"1222590007","display":"American Joint Committee on Cancer pathological N category allowable value (qualifier value)"}
tnm_pn_dict = {0:{"system":"http://snomed.info/sct","code":"1229945006","display":"pNX"},
               1:{"system":"http://snomed.info/sct","code":"1229947003","display":"pN0"},
               2:{"system":"http://snomed.info/sct","code":"1229951001","display":"pN1"},
               3:{"system":"http://snomed.info/sct","code":"1229957002","display":"pN2"}}
tnm_cn_code = {"system":"http://snomed.info/sct","code":"1222588006","display":"American Joint Committee on Cancer clinical N category allowable value (qualifier value)"}
tnm_cn_dict = {0:{"system":"http://snomed.info/sct","code":"1229966003","display":"cNX"},
               1:{"system":"http://snomed.info/sct","code":"1229967007","display":"cN0"},
               2:{"system":"http://snomed.info/sct","code":"1229973008","display":"cN1"},
               3:{"system":"http://snomed.info/sct","code":"1229978004","display":"cN2"}}
tnm_pm_code = {"system":"http://snomed.info/sct","code":"1222591006","display":"American Joint Committee on Cancer pathological M category allowable value (qualifier value)"}
tnm_pm_dict = {0:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C48740","display":"pM0 TNM Finding"},
               1:{"system":"http://snomed.info/sct","code":"1229916009","display":"pM1"},
               2:{"system":"http://snomed.info/sct","code":"1229917000","display":"pM1a"},
               3:{"system":"http://snomed.info/sct","code":"1229920008","display":"pM1b"},
               4:{"system":"http://snomed.info/sct","code":"1229923005","display":"pM1c"}}
tnm_cm_code = {"system":"http://snomed.info/sct","code":"1222587001","display":"American Joint Committee on Cancer clinical M category allowable value (qualifier value)"}
tnm_cm_dict = {0:{"system":"http://snomed.info/sct","code":"1229901006","display":"cM0"},
               1:{"system":"http://snomed.info/sct","code":"1229903009","display":"cM1"},
               2:{"system":"http://snomed.info/sct","code":"1229904003","display":"cM1a"},
               3:{"system":"http://snomed.info/sct","code":"1229907005","display":"cM1b"},
               4:{"system":"http://snomed.info/sct","code":"1229910003","display":"cM1c"}}
notrecorded = {"system":"http://snomed.info/sct", "code":"1220561009","display":"Not recorded (qualifier value)"}
yesno_dict = {0:{"system":"http://snomed.info/sct","code":"373066001","display":"Yes"},
              1:{"system":"http://snomed.info/sct","code":"260415000","display":"Not detected"},
              2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
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

#define function for histopat
def hist(hist_index,paturl,bundle,urnpat,urnint,urncase,note=np.nan):
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
    #add note
    if pd.isna(note) == False:
        hist_data["note"] = [{"text":note}]

    histopat = Condition(**hist_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(urnpat)+str(urnint)+str(urncase)).urn,"resource":histopat})

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
        "id": str("UC3-bundle-")+str(d_index["redcap_data_access_group_index"])+"-"+str(index+1),
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
                    "code":"94381002:260716003=93761005",
                    "display":"Secondary malignant neoplasm of liver where Site of primary tumour = Primary malignant neoplasm of colon"
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
    #evidence setup
    if pd.isna(d_index["diag_method_index"]) == False and d_index["diag_method_index"] in diag_method_dict:
        primary_condition_data["evidence"] = [{"code":[{"coding":[diag_method_dict[d_index["diag_method_index"]]]}]}]
    else:
        primary_condition_data["evidence"] = [{"code":[{"coding":[notrecorded]}]}]
    
    primary_condition = Condition(**primary_condition_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("condition")+str(3)).urn,"resource":primary_condition})
    
    #histopat setup
    hist(d_index["histopat_index"],paturl,bundle,d_index["record_id_index"],"histopat",3,note="Histological type of liver lesion")
    #comorb_present setup
    obs(d_index["comorb_present_index"],yesno_dict,comorb_code,paturl,bundle,d_index["record_id_index"],"comorb",3)
    #bca setup
    obs(d_index["bca_index"],yesno_dict,bca_code,paturl,bundle,d_index["record_id_index"],"bca",3)
    #biloma setup
    obs(d_index["biloma_index"],yesno_dict,biloma_code,paturl,bundle,d_index["record_id_index"],"biloma",3, note="Biloma")
    #cca setup
    obs(d_index["cca_index"],yesno_dict,cca_code,paturl,bundle,d_index["record_id_index"],"cca",3)
    #dn setup
    obs(d_index["dn_index"],yesno_dict,dn_code,paturl,bundle,d_index["record_id_index"],"dn",3)
    #fhcc setup
    obs(d_index["fhcc_index"],yesno_dict,fhcc_code,paturl,bundle,d_index["record_id_index"],"fhcc",3)
    #fhs setup
    obs(d_index["fhs_index"],yesno_dict,fhs_code,paturl,bundle,d_index["record_id_index"],"fhs",3)
    #hemangioma setup
    obs(d_index["hemangioma_index"],yesno_dict,hemangioma_code,paturl,bundle,d_index["record_id_index"],"hemangioma",3)
    #abscess setup
    obs(d_index["abscess_index"],yesno_dict,abscess_code,paturl,bundle,d_index["record_id_index"],"abscess",3)
    #adenoma setup
    obs(d_index["adenoma_index"],yesno_dict,adenoma_code,paturl,bundle,d_index["record_id_index"],"adenoma",3)
    #pha setup
    obs(d_index["pha_index"],yesno_dict,pha_code,paturl,bundle,d_index["record_id_index"],"pha",3)
    #hhcyst setup
    obs(d_index["hhcyst_index"],yesno_dict,hhcyst_code,paturl,bundle,d_index["record_id_index"],"hhcyst",3)
    #hcc setup
    obs(d_index["hcc_index"],yesno_dict,hcc_code,paturl,bundle,d_index["record_id_index"],"hcc",3)
    #cirrhosis setup
    obs(d_index["cirrhosis_index"],yesno_dict,cirrhosis_code,paturl,bundle,d_index["record_id_index"],"cirrhosis",3)
    #phl setup
    obs(d_index["phl_index"],yesno_dict,phl_code,paturl,bundle,d_index["record_id_index"],"phl",3)
    #shcyst setup
    obs(d_index["shcyst_index"],yesno_dict,shcyst_code,paturl,bundle,d_index["record_id_index"],"shcyst",3)
    #thad setup
    obs(d_index["thad_index"],yesno_dict,thad_code,paturl,bundle,d_index["record_id_index"],"thad",3,note="Presence/Absence of Transient hepatic attenuation differences (THAD)")
    #tnm_pt setup
    obs(d_index["tnm_pt_index"],tnm_pt_dict,tnm_pt_code,paturl,bundle,d_index["record_id_index"],"tnm_pt",3)
    #tnm_ct setup
    obs(d_index["tnm_ct_index"],tnm_ct_dict,tnm_ct_code,paturl,bundle,d_index["record_id_index"],"tnm_ct",3)
    #tnm_pn setup
    obs(d_index["tnm_pn_index"],tnm_pn_dict,tnm_pn_code,paturl,bundle,d_index["record_id_index"],"tnm_pn",3)
    #tnm_cn setup
    obs(d_index["tnm_cn_index"],tnm_cn_dict,tnm_cn_code,paturl,bundle,d_index["record_id_index"],"tnm_cn",3)
    #tnm_pm setup
    obs(d_index["tnm_pm_index"],tnm_pm_dict,tnm_pm_code,paturl,bundle,d_index["record_id_index"],"tnm_pm",3)
    #tnm_cm setup
    obs(d_index["tnm_cm_index"],tnm_cm_dict,tnm_cm_code,paturl,bundle,d_index["record_id_index"],"tnm_cm",3)
    
    #Making JSON readable and parsing it
    prettybundle = json.loads(bundle.json())
    parsed = json.dumps(prettybundle, indent=4)
    #Giving JSON an export name and path
    relative_path_json = "/UC3/OUTPUT/"
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