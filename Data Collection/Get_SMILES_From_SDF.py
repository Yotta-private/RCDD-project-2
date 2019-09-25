import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.six import StringIO
sio = StringIO()
w = Chem.SDWriter(sio)

#"2019-09-06-malacards_32_add_from_pubchem.sdf"
data = Chem.SDMolSupplier("2019-09-06-malacards_32_add_from_pubchem.sdf")
#len(data)
#输出数据集中分子包含的属性
#properties = data[0].GetPropNames()
#for prop in properties:
    #print(prop)
Name_smiles = [data[i].GetProp("Name") for i in range(len(data))]
smiles = [Chem.MolToSmiles(data[i]) for i in range(len(data))]

w=Chem.SDWriter("./2019-09-06-malacards_32_add_rdkit.sdf")
for i in range(len(Name_smiles)):
    m=Chem.MolFromSmiles(smiles[i])
    #print(Chem.MolToMolBlock(m)) 
    m.SetProp("_Name",Name_smiles[i])
    AllChem.Compute2DCoords(m)
    w.write(m)
w.close()

List = list(zip(smiles, Name_smiles))
datapd=pd.DataFrame(columns=["SMILES", "name"], data=List)
datapd.to_csv('./2019-09-06-malacards_32_add_rdkit.csv',encoding='utf-8')
