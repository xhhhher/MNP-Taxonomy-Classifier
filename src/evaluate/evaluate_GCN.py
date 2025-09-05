import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from matplotlib.colors import Normalize

df=pd.read_csv('../results/predictions/CMNPD2.0_test_set_predict.csv')
df["pred"] = df[["pred_0", "pred_1", "pred_2"]].idxmax(axis=1).str[-1].astype(int)
df = df.drop(columns=["pred_0", "pred_1", "pred_2"])
df_true=pd.read_csv('.../data/processed/CMNPD2.0_test_set.csv')

true_labels = df_true['labels'].values
predicted_labels = df['pred'].values

cm = confusion_matrix(true_labels, predicted_labels)
cm_percentage = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100

fig, ax = plt.subplots(figsize=(8, 6))

cax = ax.imshow(cm_percentage, interpolation='nearest', cmap='Blues', alpha=0.6)

ax.set_title('Confusion Matrix of GCN', fontsize=11,pad=15)
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
plt.savefig('.../results/figures/CMNPD2.0_test_set_ConfusionMatrix_GCN.svg',format='svg')
