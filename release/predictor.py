from __future__ import print_function
from __future__ import division

import numpy as np

#from sklearn.externals import joblib
import joblib
from sklearn import metrics
from sklearn.model_selection import train_test_split

from utils import get_fp, get_desc, normalize_desc, cross_validation_split

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
        
# class VanillaQSAR(object):
#     def __init__(self, model_instance=None, model_params=None,
#                  model_type='classifier', ensemble_size=5, normalization=False):
#         super(VanillaQSAR, self).__init__()
#         self.model_instance = model_instance
#         self.model_params = model_params
#         self.ensemble_size = ensemble_size
#         self.model = []
#         self.normalization = normalization
#         if model_type not in ['classifier', 'regressor']:
#             raise InvalidArgumentError("model type must be either"
#                                        "classifier or regressor")
#         self.model_type = model_type
#         if isinstance(self.model_instance, list):
#             assert(len(self.model_instance) == self.ensemble_size)
#             assert(isinstance(self.model_params, list))
#             assert(len(self.model_params) == self.ensemble_size)
#             for i in range(self.ensemble_size):
#                 self.model.append(self.model_instance[i](**model_params[i]))
#         else:
#             for _ in range(self.ensemble_size):
#                 self.model.append(self.model_instance(**model_params))
#         if self.normalization:
#             self.desc_mean = [0]*self.ensemble_size
#         self.metrics_type = None

#     def fit_model(self, data, cv_split='stratified'):
#         eval_metrics = []
#         x = data.x
#         if self.model_type == 'classifier' and data.binary_y is not None:
#             y = data.binary_y
#         else:
#             y = data.y
        
#         # train_x, test_x, train_y, test_y = train_test_split(x, y, test_size=0.05, random_state=9)
#         # self.model[0].fit(train_x, train_y.ravel())
#         # predicted = self.model[0].predict(test_x)
#         # if self.model_type == 'classifier':
#         #     eval_metrics.append(metrics.f1_score(test_y, predicted))
#         #     self.metrics_type = 'F1 score'
#         # elif self.model_type == 'regressor':
#         #     r2 = metrics.r2_score(test_y, predicted)
#         #     eval_metrics.append(r2)
#         #     self.metrics_type = 'R^2 score'
#         # else:
#         #     raise RuntimeError()

#         cross_val_data, cross_val_labels = cross_validation_split(x=x, y=y,
#                                                                   split=cv_split,
#                                                                   n_folds=self.ensemble_size)
#         for i in range(self.ensemble_size):
#             train_x = np.concatenate(cross_val_data[:i] +
#                                      cross_val_data[(i + 1):])
#             test_x = cross_val_data[i]
#             train_y = np.concatenate(cross_val_labels[:i] +
#                                      cross_val_labels[(i + 1):])
#             test_y = cross_val_labels[i]
#             if self.normalization:
#                 train_x, desc_mean = normalize_desc(train_x)
#                 self.desc_mean[i] = desc_mean
#                 test_x, _ = normalize_desc(test_x, desc_mean)
#             self.model[i].fit(train_x, train_y.ravel())
#             predicted = self.model[i].predict(test_x)
#             if self.model_type == 'classifier':
#                 eval_metrics.append(metrics.f1_score(test_y, predicted))
#                 self.metrics_type = 'F1 score'
#             elif self.model_type == 'regressor':
#                 r2 = metrics.r2_score(test_y, predicted)
#                 eval_metrics.append(r2)
#                 self.metrics_type = 'R^2 score'
#             else:
#                 raise RuntimeError()
#         return eval_metrics, self.metrics_type

#     def load_model(self, path):
#         # TODO: add iterable path object instead of static path 
#         self.model = []
#         for i in range(self.ensemble_size):
#             m = joblib.load(path + str(i) + '.pkl')
#             self.model.append(m)
#         if self.normalization:
#             arr = np.load(path + 'desc_mean.npy')
#             self.desc_mean = arr

#     def save_model(self, path):
#         assert self.ensemble_size == len(self.model)
#         for i in range(self.ensemble_size):
#             joblib.dump(self.model[i], path + str(i) + '.pkl')
#         if self.normalization:
#             np.save(path + 'desc_mean.npy', self.desc_mean)

#     def predict(self, objects=None, average=True, get_features=None,
#                 **kwargs):
#         objects = np.array(objects)
#         invalid_objects = []
#         processed_objects = []
#         if get_features is not None:
#             x, processed_indices, invalid_indices = get_features(objects,
#                                                                  **kwargs)
#             processed_objects = objects[processed_indices]
#             invalid_objects = objects[invalid_indices]
#         else:
#             x = objects
#         if len(x) == 0:
#             processed_objects = []
#             prediction = []
#             invalid_objects = objects
#         else:
#             prediction = []
#             for i in range(self.ensemble_size):
#                 m = self.model[i]
#                 if self.normalization:
#                     x, _ = normalize_desc(x, self.desc_mean[i])
#                 prediction.append(m.predict(x))
#             prediction = np.array(prediction)
#             if average:
#                 prediction = prediction.mean(axis=0)
#         return processed_objects, prediction, invalid_objects
