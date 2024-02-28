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
from fhir.resources.procedure import Procedure
import json

# %%
#Import data from CSV. Change formats to readable ones by dictionaries.
#absolute_path = os.path.dirname(os.getcwd())
absolute_path = os.getcwd()
relative_path_csv = "/UC4&5/CSV/UseCase45_testdata.csv"
csvpath = absolute_path+relative_path_csv
data_csv = pd.read_csv(csvpath)
data_csv = data_csv.astype({
    "sex":"Int64",
    "histopat":"Int64",
    "histopat_2":"Int64",
    "histopat_3":"Int64",
    "tnm_ct":"Int64",
    "tnm_ypt":"Int64",
    "tnm_ypn":"Int64",
    "tnm_ypm":"Int64",
    "grade":"Int64",
    "lymp_inv":"Int64",
    "perineural_inv":"Int64",
    "emvi":"Int64",
    "n_lymphnodes":"Int64",
    "n_met_ln":"Int64",
    "tumor_deposits":"Int64",
    "trg":"Int64",
    "radiotherapy":"Int64",
    "surgery":"Int64",
    "chemo":"Int64",
    "chemo_drug":"Int64",
    "neoadj":"Int64",
    "chemo_drug_2":"Int64",
    "neoadj_2":"Int64",
    "chemo_drug_3":"Int64",
    "neoadj_3":"Int64",
    "chemo_drug_4":"Int64",
    "neoadj_4":"Int64",
    "chemo_drug_5":"Int64",
    "neoadj_5":"Int64",
    "neoadj_surgery":"Int64"
})

#Convert table to dictionary containing all variables as lists
d_csv = {}
for i in range(len(data_csv.columns)):
    d_csv[data_csv.columns[i]+"_csv"] = data_csv[data_csv.columns[i]].tolist()

# %%
sex_code= {"system": "http://loinc.org","code": "76689-9","display": "Sex assigned at birth"}
sex_dict = {0: {"system":"http://loinc.org","code":"LA2-8","display":"Male"},
            1: {"system":"http://loinc.org","code":"LA3-6","display":"Female"}}
histopat_dict = {0:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8140/3","display":"Adenocarcinoma NOS"},
                 1:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8246/3","display":"Neuroendocrine Carcinoma, NOS"},
                 2:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8154/3","display":"Mixed neuroendocrine non-neuroendocrine neoplasm (MiNEN)"},
                 3:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8480/3","display":"Mucinous adenocarcinoma"},
                 4:{"system":"http://terminology.hl7.org/CodeSystem/icd-o-3","code":"8255/3","display":"Adenocarcinoma with mixed subtypes"}}
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
tnm_ypt_code = {"system":"http://snomed.info/sct","code":"1222589003","display":"American Joint Committee on Cancer pathological T category allowable value (qualifier value)"}
tnm_ypt_dict = {0:{"system":"http://snomed.info/sct","code":"1228950008","display":"pTx"},
               1:{"system":"http://snomed.info/sct","code":"1228951007","display":"pT0"},
               2:{"system":"http://snomed.info/sct","code":"1228953005","display":"pTis"},
               3:{"system":"http://snomed.info/sct","code":"1228957006","display":"pT1"},
               4:{"system":"http://snomed.info/sct","code":"1229852009","display":"pT2"},
               5:{"system":"http://snomed.info/sct","code":"1229859000","display":"pT3"},
               6:{"system":"http://snomed.info/sct","code":"1229864001","display":"pT4"},
               7:{"system":"http://snomed.info/sct","code":"1229865000","display":"pT4a"},
               8:{"system":"http://snomed.info/sct","code":"1229866004","display":"pT4b"}}
tnm_ypn_code = {"system":"http://snomed.info/sct","code":"1222590007","display":"American Joint Committee on Cancer pathological N category allowable value (qualifier value)"}
tnm_ypn_dict = {0:{"system":"http://snomed.info/sct","code":"1229945006","display":"pNX"},
               1:{"system":"http://snomed.info/sct","code":"1229947003","display":"pN0"},
               2:{"system":"http://snomed.info/sct","code":"1229951001","display":"pN1"},
               3:{"system":"http://snomed.info/sct","code":"1229957002","display":"pN2"}}
tnm_ypm_code = {"system":"http://snomed.info/sct","code":"1222591006","display":"American Joint Committee on Cancer pathological M category allowable value (qualifier value)"}
tnm_ypm_dict = {0:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C48740","display":"pM0 TNM Finding"},
               1:{"system":"http://snomed.info/sct","code":"1229916009","display":"pM1"},
               2:{"system":"http://snomed.info/sct","code":"1229917000","display":"pM1a"},
               3:{"system":"http://snomed.info/sct","code":"1229920008","display":"pM1b"},
               4:{"system":"http://snomed.info/sct","code":"1229923005","display":"pM1c"}}
grade_code = {"system":"http://snomed.info/sct","code":"258244004","display":"Tumor histopathological grade status values (tumor staging)"}
grade_dict = {0:{"system":"http://snomed.info/sct","code":"1228845001","display":"American Joint Committee on Cancer grade GX (qualifier value)"},
              1:{"system":"http://snomed.info/sct","code":"1228848004","display":"American Joint Committee on Cancer grade G1 (qualifier value)"},
              2:{"system":"http://snomed.info/sct","code":"1228850007","display":"American Joint Committee on Cancer grade G2 (qualifier value)"},
              3:{"system":"http://snomed.info/sct","code":"1228851006","display":"American Joint Committee on Cancer grade G3 (qualifier value)"}}
lymp_inv_code = {"system":"http://loinc.org","code":"LP428220-0","display":"Lymphovascular invasion"}
perineural_inv_code = {"system":"http://loinc.org","code":"92837-4","display":"Perineural invasion [Presence] in Cancer specimen"}
emvi_code = {"system":"http://loinc.org","code":"84889-5","display":"Extramural vein invasion [Presence] in Colorectal cancer specimen by Light microscopy"}
n_lymphnodes_code = {"system":"http://loinc.org","code":"21894-1","display":"Regional lymph nodes examined [#] Specimen"}
n_met_ln_code = {"system":"http://snomed.info/sct","code":"443527007","display":"Number of lymph nodes containing metastatic neoplasm in excised specimen (observable entity)"}
tumor_deposits_code = {"system":"http://snomed.info/sct","code":"630001000004109",
                       "display":"Presence of metastatic discontinuous spread of malignant neoplasm of colon to pericolic tissue (observable entity)"}
trg_code = {"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C155939","display":"Modified Ryan Scheme for Tumor Regression"}
trg_dict = {0:{"system":"http://snomed.info/sct","code":"1155705000","display":"GX: Histologic grade cannot be assessed (qualifier value)"},
            1:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C155941","display":"Tumor Regression Score 0"},
            2:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C155942","display":"Tumor Regression Score 1"},
            3:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C155943","display":"Tumor Regression Score 2"},
            4:{"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C155944","display":"Tumor Regression Score 3"}}
radiotherapy_code = {"system":"http://snomed.info/sct","code":"108290001","display":"Radiation oncology AND/OR radiotherapy (procedure)"}
surgery_code = {"system":"http://snomed.info/sct","code":"387713003","display":"Surgical procedure (procedure)"}
chemo_code = {"system":"http://snomed.info/sct","code":"400001000004103","display":"Neoadjuvant antineoplastic therapy (procedure)"}
chemo_drug_dict = {0: {"system":"http://www.nlm.nih.gov/research/umls/rxnorm","code":"194000","display":"Capecitabine"},
                   1: {"system":"http://www.nlm.nih.gov/research/umls/rxnorm","code":"4492","display":"Fluorouracil"},
                   2: {"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C63593","display":"FOLFIRI"},
                   3: {"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C11764","display":"FOLFIRINOX"},
                   4: {"system":"http://ncithesaurus-stage.nci.nih.gov","code":"C11197","display":"FOLFOX"}}
neoadj_surgery_code = {"system":"http://snomed.info/sct","code":"307474000:398201009=400001000004103,397898000=387713003",
                       "display":"Intervals of weeks and months where Start time = Neoadjuvant antineoplastic therapy, and Stop time = Surgical procedure"}

#Units and generic dictionaries
month_unit = {"value":"","unit":"month","system":"http://unitsofmeasure.org","code":"mo"}
unit_unit = {"value":"","unit":"Unit","system":"http://unitsofmeasure.org","code":"U"}
percent_unit = {"value":"","unit":"percent","system":"http://unitsofmeasure.org","code":"%"}
mgm2_unit = {"value":"","unit":"milligram per square meter","system":"http://unitsofmeasure.org","code":"mg/m2"}
mgkg_unit = {"value":"","unit":"milligram per kilogram","system":"http://unitsofmeasure.org","code":"mg/kg"}
mg_unit = {"value":"","unit":"milligram","system":"http://unitsofmeasure.org","code":"mg"}
auc2_unit = {"value":"","unit":"percent","system":"http://unitsofmeasure.org","code":""}
auc5_unit = {"value":"","unit":"percent","system":"http://unitsofmeasure.org","code":""}

#chemo_unit_dict = {0:mgm2_unit,1:mgkg_unit,2:mg_unit,3:auc2_unit,4:auc5_unit}
chemo_unit_dict = {0:mgm2_unit,1:mgkg_unit,2:mg_unit}

yesno_dict = {0:{"system":"http://snomed.info/sct","code":"373066001","display":"Yes"},
              1:{"system":"http://snomed.info/sct","code":"260415000","display":"Not detected"},
              2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
posnegun_dict = {0:{"system":"http://snomed.info/sct","code":"10828004","display":"Positive"},
                 1:{"system":"http://snomed.info/sct","code":"260385009","display":"Negative"},
                 2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown"}}
presence_dict = {0:{"system":"http://snomed.info/sct","code":"52101004","display":"Present (qualifier value)"},
                 1:{"system":"http://snomed.info/sct","code":"2667000","display":"Absent (qualifier value)"},
                 2:{"system":"http://snomed.info/sct","code":"261665006","display":"Unknown (qualifier value)"}}
identified_dict = {0:{"system":"http://loinc.org","code":"LA11902-6","display":"Not identified"},
                 1:{"system":"http://loinc.org","code":"LA9633-4","display":"Present"}}
status_dict = {0: "not-done", 1: "stopped", 2: "completed", 3: "unknown"}
status_dict_2 = {0: "not-done", 1: "completed", 2: "unknown"}
notrecorded = {"system":"http://snomed.info/sct", "code":"1220561009","display":"Not recorded (qualifier value)"}
#{"system":"","code":"","display":""}

# %%
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

#define function for procedure
def proc(proc_index,proc_code,stat_dict,paturl,bundle,urnpat,urnint,urncase,note=np.nan):
    proc_data = proc_data = {
            "code":{"coding":[proc_code]},
            "subject":{"reference":paturl}
        }
    if pd.isna(proc_index) == False and proc_index in stat_dict:
        proc_data["status"]=stat_dict[proc_index]
    else:
        proc_data["status"]="unknown"
    #add note
    if pd.isna(note) == False:
        proc_data["note"] = [{"text":note}]

    proced = Procedure(**proc_data)
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(urnpat)+str(urnint)+str(urncase)).urn,"resource":proced})

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
        "id": str("UC45-bundle-")+str(d_index["redcap_data_access_group_index"])+"-"+str(index+1),
        "type": "collection",
        "entry":[]
    }
    bundle = Bundle(**bundle_data)
    
    #Patient setup
    paturl = uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("patient")+str(45)).urn
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
    obs(d_index["sex_index"],sex_dict,sex_code,paturl,bundle,d_index["record_id_index"],"sex",45)
    
    #primary condition setup
    primary_condition_data ={
        "code":{"coding":[{
                    "system":"http://snomed.info/sct",
                    "code":"93984006",
                    "display":"Primary malignant neoplasm of rectum (disorder)"
        }]},
        "subject":{"reference":paturl},
        "bodySite":[
            {
                "coding":[
                    {
                        "system":"http://snomed.info/sct",
                        "code":"34402009",
                        "display":"Rectum structure (body structure)"
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
    bundle.entry.append({"fullUrl":uuid.uuid5(uuid.NAMESPACE_DNS,str(d_index["record_id_index"])+str("condition")+str(45)).urn,"resource":primary_condition})
    
    #histopat 1 setup
    hist(d_index["histopat_index"],paturl,bundle,d_index["record_id_index"],"histopat1",45,note="Histopathological type")
    #histopat 2 setup
    hist(d_index["histopat_2_index"],paturl,bundle,d_index["record_id_index"],"histopat2",45,note="Additional histopathological type")
    #histopat 3 setup
    hist(d_index["histopat_3_index"],paturl,bundle,d_index["record_id_index"],"histopat3",45,note="Second additional Histopathological type")
    #tnm_ct setup
    obs(d_index["tnm_ct_index"],tnm_ct_dict,tnm_ct_code,paturl,bundle,d_index["record_id_index"],"tnm_ct",45)
    #tnm_ypt setup
    obs(d_index["tnm_ypt_index"],tnm_ypt_dict,tnm_ypt_code,paturl,bundle,d_index["record_id_index"],"tnm_ypt",45)
    #tnm_ypn setup
    obs(d_index["tnm_ypn_index"],tnm_ypn_dict,tnm_ypn_code,paturl,bundle,d_index["record_id_index"],"tnm_ypn",45)
    #tnm_ypm setup
    obs(d_index["tnm_ypm_index"],tnm_ypm_dict,tnm_ypm_code,paturl,bundle,d_index["record_id_index"],"tnm_ypm",45)
    #grade setup
    obs(d_index["grade_index"],grade_dict,grade_code,paturl,bundle,d_index["record_id_index"],"grade",45)
    #lymp_inv setup
    obs(d_index["lymp_inv_index"],presence_dict,lymp_inv_code,paturl,bundle,d_index["record_id_index"],"lymp_inv",45)
    #perineural_inv setup
    obs(d_index["perineural_inv_index"],presence_dict,perineural_inv_code,paturl,bundle,d_index["record_id_index"],"perineural_inv",45)
    #emvi setup
    obs(d_index["emvi_index"],presence_dict,emvi_code,paturl,bundle,d_index["record_id_index"],"emvi",45)
    #n_lymphnodes setup
    obsquant(d_index["n_lymphnodes_index"],unit_unit,n_lymphnodes_code,paturl,bundle,d_index["record_id_index"],"n_lymphnodes",45)
    #n_met_ln setup
    obsquant(d_index["n_met_ln_index"],unit_unit,n_met_ln_code,paturl,bundle,d_index["record_id_index"],"n_met_ln",45)
    #tumor_deposits setup
    obs(d_index["tumor_deposits_index"],presence_dict,tumor_deposits_code,paturl,bundle,d_index["record_id_index"],"tumor_deposits",45)
    #trg setup
    obs(d_index["trg_index"],trg_dict,trg_code,paturl,bundle,d_index["record_id_index"],"trg",45)
    #radiotherapy setup
    proc(d_index["radiotherapy_index"],radiotherapy_code,status_dict,paturl,bundle,d_index["record_id_index"],"radiotherapy",45)
    #surgery setup
    proc(d_index["surgery_index"],surgery_code,status_dict_2,paturl,bundle,d_index["record_id_index"],"surgery",45)
    #chemo setup
    proc(d_index["chemo_index"],chemo_code,status_dict,paturl,bundle,d_index["record_id_index"],"chemo",45)
    #medicaton 1 setup
    medadm(d_index["chemo_drug_index"],chemo_drug_dict,paturl,bundle,
           d_index["record_id_index"],"med1",45,
           med_status=3,med_status_dict=status_dict,
           med_dose=d_index["chemo_dose_index"],med_unit=0,med_unit_dict=chemo_unit_dict,
           med_cycles=d_index["neoadj_index"],
           effectivetimeweeks=d_index["chemo_dur_index"],
           note="The effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates. If the effectivePeriod is the same at the start and end, no duration was provided")
    #medicaton 2 setup
    medadm(d_index["chemo_drug_2_index"],chemo_drug_dict,paturl,bundle,
           d_index["record_id_index"],"med2",45,
           med_status=3,med_status_dict=status_dict,
           med_dose=d_index["chemo_dose_2_index"],med_unit=0,med_unit_dict=chemo_unit_dict,
           med_cycles=d_index["neoadj_2_index"],
           effectivetimeweeks=d_index["chemo_dur_2_index"],
           note="The effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates. If the effectivePeriod is the same at the start and end, no duration was provided")
    #medicaton 2 setup
    medadm(d_index["chemo_drug_3_index"],chemo_drug_dict,paturl,bundle,
           d_index["record_id_index"],"med3",45,
           med_status=3,med_status_dict=status_dict,
           med_dose=d_index["chemo_dose_3_index"],med_unit=0,med_unit_dict=chemo_unit_dict,
           med_cycles=d_index["neoadj_3_index"],
           effectivetimeweeks=d_index["chemo_dur_3_index"],
           note="The effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates. If the effectivePeriod is the same at the start and end, no duration was provided")
    #medicaton 4 setup
    medadm(d_index["chemo_drug_4_index"],chemo_drug_dict,paturl,bundle,
           d_index["record_id_index"],"med4",45,
           med_status=3,med_status_dict=status_dict,
           med_dose=d_index["chemo_dose_4_index"],med_unit=0,med_unit_dict=chemo_unit_dict,
           med_cycles=d_index["neoadj_4_index"],
           effectivetimeweeks=d_index["chemo_dur_4_index"],
           note="The effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates. If the effectivePeriod is the same at the start and end, no duration was provided")
    #medicaton 5 setup
    medadm(d_index["chemo_drug_5_index"],chemo_drug_dict,paturl,bundle,
           d_index["record_id_index"],"med5",45,
           med_status=3,med_status_dict=status_dict,
           med_dose=d_index["chemo_dose_5_index"],med_unit=0,med_unit_dict=chemo_unit_dict,
           med_cycles=d_index["neoadj_5_index"],
           effectivetimeweeks=d_index["chemo_dur_5_index"],
           note="The effectivePeriod variable has anonymized dates. The duration is the difference in time between both dates. If the effectivePeriod is the same at the start and end, no duration was provided")
    #neoadj_surgery setup
    obsquant(d_index["neoadj_surgery_index"],month_unit,neoadj_surgery_code,paturl,bundle,d_index["record_id_index"],"neoadj_surgery",45)
    #Making JSON readable and parsing it
    prettybundle = json.loads(bundle.json())
    parsed = json.dumps(prettybundle, indent=4)
    #Giving JSON an export name and path
    relative_path_json = "/UC4&5/OUTPUT/"
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


