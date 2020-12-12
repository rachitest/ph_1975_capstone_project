# %%
import yaml
import json
import os

import sqlalchemy as sql
import pandas as pd
import plotly.express as px

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

def clean_data(records): # TODO REFACTOR ASAP
    ''' Using a list of dictionaries (that contains all citation data for the dataset), on a per citation basis, it extracts the following information about the citations where possible:
    title, abstract, date, authors. The extracted information is saved as a list which is then converted into a df. 
    ''' 
    for record in records:
        if record.get("PubmedArticle") != []:
            a = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"]
            if "Abstract" in record["PubmedArticle"][0]["MedlineCitation"]["Article"].keys():
                b = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            else:
                b = []
            if ("ArticleDate" in record["PubmedArticle"][0]["MedlineCitation"]["Article"].keys()) and ((record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"]) != []):
                clean_date = pd.json_normalize(record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"]).values.tolist()
                clean_date = [item for sublist in clean_date for item in sublist]
                c = "-".join(clean_date)
            else:
                c = []
            if "AuthorList" in record["PubmedArticle"][0]["MedlineCitation"]["Article"].keys():
                clean_name = pd.json_normalize(record["PubmedArticle"][0]["MedlineCitation"]["Article"]["AuthorList"])
                if "LastName" in clean_name and "ForeName" in clean_name:
                    clean_name = clean_name["LastName"] + " " + clean_name["ForeName"]
                elif "CollectiveName" in clean_name:
                    clean_name = clean_name["CollectiveName"]
                elif "ForeName" not in clean_name or "CollectiveName" not in clean_name:
                    clean_name = clean_name["LastName"]
                elif "LastName" not in clean_name or "CollectiveName" not in clean_name:
                    clean_name = clean_name["ForeName"]
                d = clean_name.values.tolist()
            else:
                d = []
            e = record["PubmedArticle"][0]["MedlineCitation"]["PMID"]
        elif record.get("PubmedArticle") == []:
            a = record["PubmedBookArticle"][0]["BookDocument"]["Book"]["BookTitle"]
            if "Abstract" in record["PubmedBookArticle"][0]["BookDocument"].keys():
                b = record["PubmedBookArticle"][0]["BookDocument"]["Abstract"]["AbstractText"]
            else:
                b = []
            if "PubDate" in record["PubmedBookArticle"][0]["BookDocument"]["Book"].keys():
                clean_date = pd.json_normalize(record["PubmedBookArticle"][0]["BookDocument"]["Book"]["PubDate"]).values.tolist()
                clean_date = [item for sublist in clean_date for item in sublist]
                c = "-".join(clean_date)
            else:
                c = []
            if ("AuthorList" in record["PubmedBookArticle"][0]["BookDocument"]["Book"].keys()) and (record["PubmedBookArticle"][0]["BookDocument"]["AuthorList"] != []):
                clean = record["PubmedBookArticle"][0]["BookDocument"]["AuthorList"]
                clean_name = [item for sublist in clean for item in sublist]
                clean_name = pd.json_normalize(clean_name)
                if "LastName" in clean_name and "ForeName" in clean_name:
                    clean_name = clean_name["LastName"] + " " + clean_name["ForeName"]
                elif "CollectiveName" in clean_name:
                    clean_name = clean_name["CollectiveName"]
                elif "ForeName" not in clean_name or "CollectiveName" not in clean_name:
                    clean_name = clean_name["LastName"]
                elif "LastName" not in clean_name or "CollectiveName" not in clean_name:
                    clean_name = clean_name["ForeName"]
                d = clean_name.values.tolist()
            else:
                d = []
            e = record["PubmedBookArticle"][0]["BookDocument"]["PMID"]

        v = [e,a,b,c,d]

        data_tmp = pd.DataFrame(v).transpose().rename(columns={0:"pmid",1:"title",2:"abstract",3:"date",4:"author(s)"})

        if records.index(record) == 0:
            data = data_tmp
        else:
            data = pd.concat([data,data_tmp])

    return data

def sqlite_out(clean_records):
    engine = sql.create_engine('sqlite:///HIV_Records.db', echo=False)
    clean_records.to_sql("HIV_Records", con=engine, if_exists='replace', index=False)

def sql_author_query(author_name):
    engine = sql.create_engine('sqlite:///HIV_Records.db', echo=False)
    sql_df = pd.read_sql(f"select * from HIV_Records as h where h.'author(s)' like '%{author_name}%'", con=engine)
    
    return sql_df

def keep_cleaning(df):
    df = df.reset_index(drop=True)
    df.pmid = df.pmid.astype(int)
    df.date = pd.to_datetime(df.date, format='%Y-%m-%d', utc=False, errors="coerce").dt.date()

    columns = ["title", "abstract"]
    for column in columns:
        df[column] = [','.join(item) if isinstance(item, list) else item for item in df[column]]

    return df

def csv_bnb(file_path):
    well_rested_csv = pd.read_csv(file_path, dtype={'pmid': int}, parse_dates=['date'])
    well_rested_csv["date"] = pd.to_datetime(well_rested_csv.date, format='%Y-%m-%d', utc=False).dt.date
    change = ["'","[","]", "nan", ""]
    for i in range(len(well_rested_csv)):
        for j in change:
            well_rested_csv["author(s)"][i] = well_rested_csv["author(s)"][i].replace(j, '')
    
    return well_rested_csv

def draw_graph(df, graph_type='line'):
    df_copy = df.copy()
    df_copy["month"] = pd.to_datetime(df_copy["date"]).dt.month_name()
    df_copy["day"] = pd.to_datetime(df_copy["date"]).dt.day

    months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    month = pd.DataFrame(df_copy["month"].value_counts())
    month.index = pd.CategoricalIndex(month.index, categories=months, ordered=True)
    month = month.sort_index()

    if graph_type.lower() == "line":
        month_fig = px.line(month, x=month.index, y="month", 
                            labels= {
                                "index":"Months",
                                "month":"Number of Publications"
                            },
                            title = "Monthly Trend for HIV Publications"
                        )
    elif graph_type.lower() == "bar":
        month_fig = px.bar(month, x=month.index, y="month", 
                            labels= {
                                "index":"Months",
                                "month":"Number of Publications"
                            },
                            title = "Monthly Trend for HIV Publications"
                        )
    elif graph_type.lower() == "both":
        month_fig = px.line(month, x=month.index, y="month", 
                            labels= {
                                "index":"Months",
                                "month":"Number of Publications"
                            },
                            title = "Monthly Trend for HIV Publications"
                        )
        month_fig.add_bar(x=month.index, y=month["month"])
        month_fig.layout.update(showlegend=False) 

    return month_fig.show()

def summary_stats(df, calendar_month):
    df_copy = df.copy()
    df_copy["month"] = pd.to_datetime(df_copy["date"]).dt.month_name()
    df_copy["day"] = pd.to_datetime(df_copy["date"]).dt.day
    stats_df = df_copy[df_copy["month"] == calendar_month.title()]["day"].value_counts().describe().to_frame()
    stats_df = stats_df.drop("count")

    return stats_df

# %%
if __name__ == '__main__':
    with open("apikeys.yaml", "r") as yamlfile: # YO DON'T RUN THIS RN
        keys = yaml.load(yamlfile, Loader=yaml.FullLoader)
        print("Read Successful")
    
    email = "rachit.sabharwal@uth.tmc.edu"
    search = "HIV"
    
    hiv_pmids = get_pmid(contact=email, key=keys["apikeys"]["ncbikey"]["key"], term=search, mindate="2020/01/01", maxdate="2020/09/01") # YO DON"T RUN THIS RN
    hiv_records = get_data(pmid_list=hiv_pmids, contact=email, key=keys["apikeys"]["ncbikey"]["key"]) # YO DON'T RUN THIS RN
    
    with open('hiv_records.json', 'w') as outfile: # YO DON'T RUN THIS RN
        json.dump(hiv_records, outfile)

    with open('D:\Dell_Desktop\Documents\Python Projects\ph_1975_capstone_project\webapp\hiv_records.json', 'r') as outfile:
        hiv_records = json.load(outfile)
    
    hiv_clean = clean_data(hiv_records)
    hiv_clean = keep_cleaning(hiv_clean)

    if not os.path.exists("hiv_records_clean.csv"):
        hiv_clean.to_csv("hiv_records_clean.csv", index=False, date_format='%Y-%m-%d')
    elif os.path.exists("hiv_records_clean.csv"):
        print("Your CSV is already up to date!")

    hiv_csv = csv_bnb("hiv_records_clean.csv")
    sqlite_out(hiv_csv)
    sql_df = sql_author_query("Julie")
    print(sql_df.head())

    draw_graph(hiv_csv)
    january_summary = summary_stats(hiv_csv, "january")