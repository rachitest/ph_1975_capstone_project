# %%
import yaml
import time

import pandas as pd

from Bio import Entrez

# %%
# Handles the eSearch endpoint for Entrez
def get_pmid(contact, key, term, **dates):
    ''' Using the Entrez search term, it queries the eSearch endpoint of the Entrez api to retrieve the corresponding pmids and join them to the input df. '''
    for_efetch = []
    Entrez.email = contact
    Entrez.api_key = key
    
    # Get total number of records
    handle = Entrez.esearch(db='pubmed', term=term, retmax=1, mindate=dates.get("mindate"), maxdate=dates.get("maxdate"))
    record = Entrez.read(handle)
    count = int(record['Count'])

    # Get all pmids with updated retmax
    handle = Entrez.esearch(db='pubmed', term=term, retmax=count, mindate=dates.get("mindate"), maxdate=dates.get("maxdate"))
    record = Entrez.read(handle)
    for_efetch.append(record['IdList'])

    # Change output from being a 1 item list
    for_efetch = pd.Series(for_efetch[0]).str.split(pat=",", expand=True).values.tolist()

    return for_efetch

# Handles the eFetch endpoint for Entrez
def get_data(pmid_list, contact, key):
    ''' Using the pmids, it queries the eFetch endpoint to retrieve the details for the corresponding citation as a list of dictionaries. ''' 
    to_clean = []
    counter = 0
    for i in range(len(pmid_list)):
            Entrez.email = contact
            Entrez.api_key = key
            handle = Entrez.efetch(db='pubmed', id=pmid_list[i], retmode='xml')
            record = Entrez.read(handle)
            to_clean.append(record)
            if counter == 600:
                print(f"Number of records retrieved is {len(to_clean)}")
                time.sleep(60)
                counter = 0
            
            counter += 1

    return to_clean

def clean_data(records):
    ''' Using a list of dictionaries (that contains all citation data for the dataset), on a per citation basis, it extracts the following information about the citations where possible:
    title, abstract, mesh terms, referenced chemicals. The extracted information is saved as a list which is then converted into a df. 
    ''' 
    for record in records:
        if record.get("PubmedArticle") != []:
            a = record['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
            if 'Abstract' in record['PubmedArticle'][0]['MedlineCitation']['Article'].keys():
                b = record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
            else:
                b = []
            if 'MeshHeadingList' in record['PubmedArticle'][0]['MedlineCitation'].keys():
                c = record['PubmedArticle'][0]['MedlineCitation']['MeshHeadingList']
            else:
                c = []
            if 'ChemicalList' in record['PubmedArticle'][0]['MedlineCitation'].keys():
                d = record['PubmedArticle'][0]['MedlineCitation']['ChemicalList']
            else:
                d = []
            e = record['PubmedArticle'][0]['MedlineCitation']['PMID']
        elif record.get("PubmedArticle") == []:
            a = record['PubmedBookArticle'][0]['BookDocument']['Book']['BookTitle']
            b = record['PubmedBookArticle'][0]['BookDocument']['Abstract']['AbstractText']
            c = record['PubmedBookArticle'][0]['BookDocument']['Sections']
            d = []
            e = record['PubmedBookArticle'][0]['BookDocument']['PMID']

        v = [e,a,b,c,d]

        data_tmp = pd.DataFrame(v).transpose().rename(columns={0:'pmid',1:'title',2:'abstract',3:'mesh',4:'chemicals'})

        if records.index(record) == 0:
            data = data_tmp
        else:
            data = pd.concat([data,data_tmp])

    return data

# %%
if __name__ == '__main__':
    with open("apikeys.yaml", "r") as yamlfile:
        keys = yaml.load(yamlfile, Loader=yaml.FullLoader)
        print("Read Successful")
    email = "rachit.sabharwal@uth.tmc.edu"
    search = "HIV"
    hiv_pmids = get_pmid(contact=email, key=keys["apikeys"]["ncbikey"]["key"], term=search, mindate="2020/01/01", maxdate="2020/09/01")
    hiv_records = get_data(pmid_list=hiv_pmids, contact=email, key=keys["apikeys"]["ncbikey"]["key"])
# %%
