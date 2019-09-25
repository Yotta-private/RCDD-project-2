#  CAS number transforme into SDF files using web information (batch)
from __future__ import print_function
import urllib
import re
import string
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
#from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.six import StringIO
sio = StringIO()
w = Chem.SDWriter(sio)

def GetMolFromCAS(casid=""):
    """
    Downloading the molecules from http://www.chemnet.com/cas/ by CAS ID (casid).
    CAS ID transforme into Inchi format
    then, using rdkit change Inchi number into SMILES
    """
    casid=str.strip(casid)
    localfile=urllib.request.urlopen('http://www.chemnet.com/cas/supplier.cgi?terms='+casid+'&l=&exact=dict')
    temp=localfile.readlines()
    html = []
    for t in temp:
        t_utf = t.decode('utf-8')# 统一编码方式
        html.append(t_utf)
    # b'    <td align="left">InChI=1/C4H13NO7P2/c5-3-1-2-4(6,13(7,8)9)14(10,11)12/h6H,1-3,5H2,(H2,7,8,9)(H2,10,11,12)</td>\r\n'
    for i in html:
        if re.findall('InChI=',i)==['InChI=']:
            k=i.split('    <td align="left">')
            kk=k[1].split('</td>\r\n')
            if kk[0][0:5]=="InChI":
                res=kk[0]    
            else:
                res="None"
    localfile.close()
    mol = Chem.inchi.MolFromInchi(res)
    smiles = Chem.MolToSmiles(mol)
    #mol=pybel.readstring('inchi',string.strip(res))
    #smile=mol.write('smi')
   # return string.strip(smile)
   # return smiles
    return smiles

Data = pd.read_csv("MalaCards_333_records_with_cas_20190906.csv")
Name = list(Data["Name"])
CAS = list(Data["CAS Number"])

SMILES = []
Error = [] # 记录要删除的名字
# 300 个分批做
length = len(CAS)
for i in range(0, 51):
    try:
        casid= CAS[i]
        smiles = GetMolFromCAS(casid)
        SMILES.append(smiles)
    except:
        print ("the number is invalid!!!, skip!!!")
        print(casid)
        Error.append(i)

for i in range(51, 101):
    try:
        casid= CAS[i]
        smiles = GetMolFromCAS(casid)
        SMILES.append(smiles)
    except:
        print ("the number is invalid!!!, skip!!!")
        print(casid)
        Error.append(i)

for i in range(101, 151):
    try:
        casid= CAS[i]
        smiles = GetMolFromCAS(casid)
        SMILES.append(smiles)
    except:
        print ("the number is invalid!!!, skip!!!")
        print(casid)
        Error.append(i)

for i in range(151, 201):
    try:
        casid= CAS[i]
        smiles = GetMolFromCAS(casid)
        SMILES.append(smiles)
    except:
        print ("the number is invalid!!!, skip!!!")
        print(casid)
        Error.append(i)

for i in range(201, 251):
    try:
        casid= CAS[i]
        smiles = GetMolFromCAS(casid)
        SMILES.append(smiles)
    except:
        print ("the number is invalid!!!, skip!!!")
        print(casid)
        Error.append(i)

for i in range(251, 301):
    try:
        casid= CAS[i]
        smiles = GetMolFromCAS(casid)
        SMILES.append(smiles)
    except:
        print ("the number is invalid!!!, skip!!!")
        print(casid)
        Error.append(i)

for i in range(301, length):
    try:
        casid= CAS[i]
        smiles = GetMolFromCAS(casid)
        SMILES.append(smiles)
    except:
        print ("the number is invalid!!!, skip!!!")
        print(casid)
        Error.append(i)

#len(SMILES)
#len(Error)
#length
# 总共333个分子，爬到了268个，还有65个各种报错
smiles_to_name = Name.copy()
#len(smiles_to_name) # 333
smiles_to_name_copy =  Name.copy()
for index in Error:
    element = smiles_to_name_copy[index]
    smiles_to_name.remove(element)
len(smiles_to_name)

w=Chem.SDWriter("./malacards_268_with_cas_20190906.sdf")
for i in range(len(smiles_to_name)):
    m=Chem.MolFromSmiles(SMILES[i])
    #print(Chem.MolToMolBlock(m)) 
    m.SetProp("_Name",smiles_to_name[i])
    AllChem.Compute2DCoords(m)
    w.write(m)
w.close()

List = list(zip(SMILES, smiles_to_name))
datapd=pd.DataFrame(columns=["SMILES", "name"], data=List)
datapd.to_csv('./malacards_268_with_cas_20190906.csv',encoding='utf-8')
