import pandas as pd

from Bio import Entrez

def get_email():
    while True:
        contact = str(input("Enter NCBI contact email address (required): "))
        try:
            if '@uth.tmc.edu' not in contact:
                raise ValueError('Email address is required!')
        except ValueError as e:
            print(e)
            continue
        else:
            break

    return contact

def get_apikey():   
    key = str(input("Enter your NCBI api key to increase rate: "))
    if key == '':
        print("You will be rate limited to 3 calls per second!")
    else:
        print("You will be rate limited to 10 calls per second!")

    return key

# Handles the eSearch endpoint for Entrez
def get_pmid(df, contact, key):
    ''' Using the Entrez search term, it queries the eSearch endpoint of the Entrez api to retrieve the corresponding pmids and join them to the input df. '''
    v = []
    for i in range(len(df)):
        Entrez.email = contact
        Entrez.api_key = key
        handle = Entrez.esearch(db='pubmed', term= df.entrezTerm.iloc[i], retmax=1)
        record = Entrez.read(handle)
        v.append(record['IdList'])

    v = pd.DataFrame(v, columns=['pmid'])
    v['pmid'] = v['pmid'].apply(str).replace(r'\b[0-9]{3}\b', 'None', regex=True)
    
    df = df.join(v, how='left')

    return df


# Handles the eSummary endpoint for Entrez
def get_doi(df, contact, key):
    ''' Using the pmids, it queries the eSummary endpoint to retrieve the corresponding dois and join them to the input df. ''' 
    v = []
    for i in range(len(df)):
        if not df.pmid.iloc[i] == 'None':
            Entrez.email = contact
            Entrez.api_key = key
            handle = Entrez.esummary(db='pubmed', id= df.pmid.iloc[i], retmax=1)
            record = Entrez.read(handle)
            info = record[0]['ArticleIds']
            v.append(info)

    x = pd.DataFrame(v)
    x = x.rename(columns={'pubmed':'pmid'})
    x.pmid = x.pmid.apply(str)
    x.doi = x.doi.apply(str)
    x.pmid = x.pmid.str.replace(r'([^0-9\s]+?)','', regex=True)

    a = x.set_index('pmid')
    df = df.set_index('pmid')
    df = df.join(a.doi)
    df = df.drop_duplicates()
    df['pmid'] = df.index
    df = df.set_index('index')

    return df


# Handles the eFetch endpoint for Entrez
def get_data(df, contact, key):
    ''' Using the pmids, it queries the eFetch endpoint to retrieve the details for the corresponding citation as a list of dictionaries. ''' 
    v = []
    for i in range(len(df)):
        if not df.pmid.iloc[i] == 'None':
            Entrez.email = contact
            Entrez.api_key = key
            handle = Entrez.efetch(db='pubmed', id=df.pmid.iloc[i], retmode='xml')
            record = Entrez.read(handle)
            v.append(record)

    return v
