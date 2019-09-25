import re 
import os
import string
import pandas as pd
import numpy as np

# <a href="/doi/10.1021/acs.jafc.7b03742">
temp = open("acs_Osteoporosis.txt", encoding='UTF-8')
txt = temp.readlines()

pattern = '<a href="/doi/(.*?)">'
DOI = []
for i in txt:
    p = re.compile(pattern)
    results = p.findall(i)
   # if len(results) != 0:
    DOI = DOI +results
len(DOI)

# len(DOI) 600
datapd=pd.DataFrame(data=DOI)
datapd.to_csv('./acs_Osteoporosis_20190911.csv',encoding='utf-8')
