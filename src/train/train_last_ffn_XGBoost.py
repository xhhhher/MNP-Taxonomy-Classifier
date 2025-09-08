import pandas as pd
import pickle
import xgboost as xgb
from sklearn.metrics import confusion_matrix, f1_score
from sklearn.model_selection import train_test_split

fp = pd.read_csv('last_FFN.csv')
data = pd.read_csv('../../data/processed/data_cmnpd_after2000.csv')

X = fp.iloc[:, :] 
y = data['labels'] 

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

model_xgb = xgb.XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
model_xgb.fit(X_train, y_train)

pickle.dump(model_xgb, open('model_XGB_after2000_all.pkl', 'wb'), protocol=4)