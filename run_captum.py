import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
import pandas as pd
import numpy as np
from torch.nn.utils.rnn import pad_sequence
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt
from captum.attr import IntegratedGradients

class GeneDataset(Dataset):
    def __init__(self, filepath, gene_subset=None, sep="\t"):
        self.data = pd.read_csv(filepath, sep=sep)
        self.data["h3k4me3"] = (self.data["h3k4me3"] - self.data["h3k4me3"].min()) / (self.data["h3k4me3"].max() - self.data["h3k4me3"].min())
        self.genes = self.data.groupby("gene_id")
        self.gene_ids = [gene for gene, group in self.genes if len(group) >= 50]
        if gene_subset:
            self.gene_ids = [g for g in self.gene_ids if g in gene_subset]

    def __len__(self):
        return len(self.gene_ids)

    def __getitem__(self, idx):
        gene_id = self.gene_ids[idx]
        gene_data = self.genes.get_group(gene_id)
        start_pos = gene_data["start_pos"].iloc[0]
        bins = gene_data[["bin_start", "bin_end", "h3k4me3"]].copy()
        TSS_bin_start = round((start_pos + 800) / 200) * 200
        promoter_bins = bins[(bins["bin_start"] >= TSS_bin_start - 1000) & (bins["bin_end"] < TSS_bin_start)]
        promoter_label = promoter_bins["h3k4me3"].values
        promoter_label = np.pad(promoter_label, (max(0, 5 - len(promoter_label)), 0), 'constant', constant_values=0)
        promoter_label = promoter_label[-5:]
        promoter_label = torch.tensor(promoter_label, dtype=torch.float32)
        gene_body_bins = bins[bins["bin_start"] >= TSS_bin_start]
        gene_body_values = torch.tensor(gene_body_bins["h3k4me3"].values, dtype=torch.float32)
        return gene_body_values, promoter_label

def collate_fn(batch):
    sequences, labels = zip(*batch)
    padded_sequences = pad_sequence(sequences, batch_first=True, padding_value=0)
    labels = torch.stack(labels)
    return padded_sequences, labels

#Model
class H3K4me3CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv1d(1, 16, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(16, 32, kernel_size=3, padding=1)
        self.global_pool = nn.AdaptiveAvgPool1d(1)
        self.fc1 = nn.Linear(32, 16)
        self.fc2 = nn.Linear(16, 5)
        self.dropout = nn.Dropout(0.4)
        self.relu = nn.ReLU()

    def forward(self, x):
        x = x.unsqueeze(1)  # [B, 1, L]
        x = self.relu(self.conv1(x))
        x = self.dropout(x)
        x = self.relu(self.conv2(x))
        x = self.dropout(x)
        x = self.global_pool(x)  # [B, 32, 1]
        x = torch.flatten(x, 1)  # [B, 32]
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return torch.sigmoid(x)  # [B, 5]

device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
filepath = "morethan50_nopromoter_nogenebody_cond_34yr_h3k4me3_heartLV.csv"
batch_size = 64

# Load fold 1 validation set
full_dataset = GeneDataset(filepath, sep="\t")
all_gene_ids = full_dataset.gene_ids
kf = KFold(n_splits=5, shuffle=True, random_state=42)
splits = list(kf.split(all_gene_ids))
_, val_idx = splits[0]
val_ids = [all_gene_ids[i] for i in val_idx]
val_dataset = GeneDataset(filepath, gene_subset=val_ids)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, collate_fn=collate_fn)

# Find max sequence length
max_len = max(inputs.shape[1] for inputs, _ in val_loader)

# Load model
model = H3K4me3CNN().to(device)
model.load_state_dict(torch.load("kfold300epochs_34yr_fold4.pt", map_location=device))
model.eval()

# Attribution method
ig = IntegratedGradients(model)

# Store all promoter bin attributions
all_avg_attributions = []

# Loop over promoter bins
for target_promoter_bin in range(5):
    print(f"Generating attribution for promoter bin {target_promoter_bin + 1}...")

    all_attributions = []

    for inputs, _ in val_loader:
        inputs = inputs.to(device)
        baseline = torch.zeros_like(inputs)
        attributions, _ = ig.attribute(inputs, baselines=baseline, target=target_promoter_bin, return_convergence_delta=True)
        pad_len = max_len - attributions.shape[1]
        if pad_len > 0:
            attributions = torch.nn.functional.pad(attributions, (0, pad_len), value=0)
        all_attributions.append(attributions.detach().cpu())

    all_attributions = torch.cat(all_attributions, dim=0)  # [N, max_len]
    avg_attr = all_attributions.mean(dim=0).numpy()  # [max_len]
    all_avg_attributions.append(avg_attr)

##plotting
plt.figure(figsize=(12, 6))
colors = ['orange', 'green', 'red', 'purple']  # For bins 2-5

for i in range(1, 5):  # Skip bin 1 (index 0)
    avg_attr = all_avg_attributions[i]
    plt.plot(np.arange(len(avg_attr)), avg_attr, label=f'Promoter Bin {i+1}', color=colors[i - 1])

plt.title('Integrated Gradients Attribution', fontsize=18)
plt.xlabel('Gene Body Bin', fontsize=16)
plt.ylabel('Attribution Score', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.xlim(0, 199)
plt.savefig('ig_h3k4me3_34yr_largefont.png', dpi=300)
plt.close()


