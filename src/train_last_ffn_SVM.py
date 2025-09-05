import pandas as pd
import pickle
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix, f1_score
from sklearn.model_selection import train_test_split

fp = pd.read_csv('last_FFN.csv')
data = pd.read_csv('../data/processed/data_cmnpd_after2000.csv')

X = fp.iloc[:, :] 
y = data['labels']  

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

clf = SVC(kernel='linear', probability=True, C=0.01)
clf.fit(X_train, y_train)

pickle.dump(clf, open('model_SVM_last_FFN.pkl', 'wb'), protocol=4)