import pandas as pd
import numpy as np
import os
import subprocess

data_path = "../../data/processed/data_cmnpd_after2000.csv"
output_dir = "../../results/models/random_split_10"   #random_split_3 random_split_5

n_repeats = 10  #3,5
total_samples = 18402
random_seed = 42

data = pd.read_csv(data_path)

label0_data = data[data["labels"] == 0].copy()
label1_data = data[data["labels"] == 1].copy()
label2_data = data[data["labels"] == 2].copy()

for repeat in range(n_repeats):
    print(f" {repeat + 1}/{n_repeats} repeats")

    # ========== label=0: 20%train+val,80%test ==========
    label0_shuffled = label0_data.sample(frac=1, random_state=random_seed + repeat)
    n_label0 = len(label0_shuffled)
    n_trainval_0 = int(n_label0 * 0.2)

    label0_trainval = label0_shuffled.iloc[:n_trainval_0].copy()
    label0_test = label0_shuffled.iloc[n_trainval_0:].copy()

    # 10% train+val --> val
    n_val_0 = int(len(label0_trainval) * 0.1)
    label0_val = label0_trainval.sample(n=n_val_0, random_state=repeat + 123)
    label0_train = label0_trainval.drop(label0_val.index)

    # ========== label=1,2: 100% train + val ==========
    label1_val = label1_data.sample(frac=0.1, random_state=repeat + 222)
    label2_val = label2_data.sample(frac=0.1, random_state=repeat + 333)

    label1_train = label1_data.drop(label1_val.index)
    label2_train = label2_data.drop(label2_val.index)

    label0_train["splits"] = "train"
    label0_val["splits"] = "val"
    label0_test["splits"] = "test"

    label1_train["splits"] = "train"
    label1_val["splits"] = "val"
    label2_train["splits"] = "train"
    label2_val["splits"] = "val"

    full_data = pd.concat([
        label0_train, label0_val, label0_test,
        label1_train, label1_val,
        label2_train, label2_val
    ], axis=0)

    split_file = os.path.join(output_dir, f"chemprop_data_random_repeat{repeat}.csv")
    model_dir = os.path.join(output_dir, f"model_repeat{repeat}")
    os.makedirs(model_dir, exist_ok=True)

    full_data.to_csv(split_file, index=False)

    train_command = [
        'chemprop', 'train',
        '--data-path', split_file,
        '--splits-column', 'splits',
        '--target-columns', 'labels',
        '--task-type', 'multiclass',
        '--multiclass-num-classes', '3',
        '--epochs', '200',
        '--batch-size', '32',
        '--output-dir', model_dir,
        '--num-workers', '0',
        '--ffn-hidden-dim', '1100',
        '--ffn-num-layers', '3',
        '--message-hidden-dim', '1100',
        '--depth', '3'
    ]

    subprocess.run(train_command)