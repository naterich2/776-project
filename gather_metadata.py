import requests
import re
import pandas as pd
samples = pd.read_csv('sample_numbers.csv')
records = []
for sample in samples.iterrows():
    data = requests.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+sample[1]['id'])
    record = {}
    for line in data.iter_lines():
        if re.search('diagnosis',str(line)):
            match = str(line)
            match = match.replace('<td style="text-align: justify">','')
            match = match.replace('</td>','')
            match = match.replace('<br>',',')
            match = match.replace("b'",'')
            match = match.replace("'",'')
            for attribute in match.split(',')[:-1]:
                key,value = attribute.split(': ')
                record[key] = value
            break
    record['id'] = sample[1]['id']
    records.append(record)
dataframe = pd.DataFrame.from_dict(records)
dataframe['name'] = samples['name']
dataframe['diagnosis'] = dataframe['diagnosis'].str.replace(r'([A-Za-z\- ]+)\(([A-Z\-]{3,5})\)',r'\2',regex=True)
dataframe.loc[:,'tissue'] = 'BA9'
dataframe.to_csv('samples.csv')

