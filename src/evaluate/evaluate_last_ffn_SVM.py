import pandas as pd
import pickle
from sklearn.metrics import confusion_matrix, f1_score
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

fp_path = "CMNPD2.0_test_set_last_FFN.csv"
mol_path = ".../data/processed/CMNPD2.0_test_set.csv"

fp = pd.read_csv(fp_path)
data = pd.read_csv(mol_path)

X = fp.iloc[:, :]
y = data['labels']

model_path = ".../results/models/model_SVM_last_FFN.pkl"

with open(model_path, 'rb') as f:
    clf = pickle.load(f)
pred_SVM = clf.predict(X)
matrix_SVM= confusion_matrix(y, pred_SVM)


accuracy_SVM = matrix_SVM.diagonal() / matrix_SVM.sum(axis=1)

f1_SVM = f1_score(y, pred_SVM, average=None)

print("SVM Accuracy (per class) on test set:", accuracy_SVM)
print("SVM F1 Score (per class) on test set:", f1_SVM)

results = pd.DataFrame({
    'True Label': y,
    'SVM Prediction': pred_SVM
})


true_labels = results['True Label'].values
predicted_labels = results['SVM Prediction'].values

cm = confusion_matrix(true_labels, predicted_labels)

cm_percentage = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100

fig, ax = plt.subplots(figsize=(8, 6))

cax = ax.imshow(cm_percentage, interpolation='nearest', cmap='Blues', alpha=0.6)

ax.set_title('Confusion Matrix of last_FFN+SVM', fontsize=11,pad=15)
ax.set_xlabel('Predicted Label', labelpad=20)
ax.set_ylabel('True Label', labelpad=20)  

ax.set_xticks(np.arange(3))
ax.set_yticks(np.arange(3))
ax.set_xticklabels(['Animalia', 'Bacteria', 'Fungi'])
ax.set_yticklabels(['Animalia', 'Bacteria', 'Fungi'])

for i in range(3):
    for j in range(3):
        ax.text(j, i, f'{cm_percentage[i, j]:.2f}%', ha='center', va='center', color='black', fontsize=10)

true_totals = cm.sum(axis=1)  
for i in range(3):
    ax.text(-0.78, i+0.2, f'(Total: {true_totals[i]})', ha='center', va='center', color='black', fontsize=9)

cbar = fig.colorbar(cax, ticks=[0, 25, 50, 75, 100])
cbar.ax.set_yticklabels([f'{x:.0f}%' for x in cbar.get_ticks()])

plt.tight_layout()
plt.savefig('.../results/figures/CMNPD2.0_test_set_ConfusionMatrix_SVM.svg',format='svg')