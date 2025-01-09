from __future__ import print_function
from __future__ import division

import numpy as np

#from sklearn.externals import joblib
import joblib
from sklearn import metrics
from sklearn.model_selection import train_test_split

from utils import get_fp, get_desc, normalize_desc, cross_validation_split

from rdkit import Chem
from rdkit.Chem import AllChem

def get_Morgan_fingerprint_frequency(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        info={}
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,3,nBits=2048,bitInfo=info)
        key=[key for key in info.keys()]
        frq=[len(value) for value in info.values()]
        new_info={key: value for key, value in zip(key, frq)}
        mff=[0] * len(list(fp))
        for key, value in new_info.items():
            mff[key] = value
    except:
        mff=None
    return mff

class Predictor():
    def __init__(self,model):
        self.model=model

    def get_Morgan_fingerprint_frequency(self,smi):
        try:
            mol = Chem.MolFromSmiles(smi)
            info={}
            fp = AllChem.GetMorganFingerprintAsBitVect(mol,3,nBits=2048,bitInfo=info)
            key=[key for key in info.keys()]
            frq=[len(value) for value in info.values()]
            new_info={key: value for key, value in zip(key, frq)}
            mff=[0] * len(list(fp))
            for key, value in new_info.items():
                mff[key] = value
        except:
            mff=None
        return mff 
    
    def get_predict_for_new_polymer(self,smi):
        fp=self.get_Morgan_fingerprint_frequency(smi)
        fp=np.array(fp)
        result=rf.predict(fp.reshape(1, -1))
        return smi,result
    
    def get_predict_for_new_polymers(self,smi):
        fp_list=[]
        smi_list=[]
        nan_smiles=[]
        for i in smi:
            if Chem.MolFromSmiles(i):
                fp=self.get_Morgan_fingerprint_frequency(i)
                #print(fp)
                fp_list.append(fp) 
                smi_list.append(i)
            else:
                nan_smiles.append(i)
        if len(fp_list)>0:
            result=rf.predict(np.array(fp_list))
        else:
            result=-1
        return smi_list,result,nan_smiles
        
