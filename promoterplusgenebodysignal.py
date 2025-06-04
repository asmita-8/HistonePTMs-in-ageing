import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from torch.utils.data import Dataset, DataLoader
from torch.nn.utils.rnn import pad_sequence
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt

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
        promoter_signal = promoter_bins["h3k4me3"].values
        promoter_signal = np.pad(promoter_signal, (max(0, 5 - len(promoter_signal)), 0), 'constant', constant_values=0)
        promoter_signal = promoter_signal[-5:]
        promoter_signal = torch.tensor(promoter_signal, dtype=torch.float32)

        gene_body_bins = bins[bins["bin_start"] >= TSS_bin_start]
        gene_body_signal = torch.tensor(gene_body_bins["h3k4me3"].values, dtype=torch.float32)

        return gene_body_signal, promoter_signal, gene_id

def collate_fn(batch):
    gene_bodies, promoters, gene_ids = zip(*batch)
    gene_bodies_padded = pad_sequence(gene_bodies, batch_first=True, padding_value=0)
    promoters_stacked = torch.stack(promoters)
    return gene_bodies_padded, promoters_stacked, gene_ids

# Model
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
        x = x.unsqueeze(1)
        x = self.relu(self.conv1(x))
        x = self.dropout(x)
        x = self.relu(self.conv2(x))
        x = self.dropout(x)
        x = self.global_pool(x)
        x = torch.flatten(x, 1)
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return torch.sigmoid(x)

# prediction and averaging Function 
def predict_and_average(filepath, model_prefix, device):
    dataset = GeneDataset(filepath)
    all_gene_ids = dataset.gene_ids
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    gene_pred_dict = {}
    gene_body_dict = {}

    for fold, (train_idx, _) in enumerate(kf.split(all_gene_ids)):
        train_ids = [all_gene_ids[i] for i in train_idx]
        train_dataset = GeneDataset(filepath, gene_subset=train_ids)
        train_loader = DataLoader(train_dataset, batch_size=64, shuffle=False, collate_fn=collate_fn)

        model = H3K4me3CNN().to(device)
        model.load_state_dict(torch.load(f"{model_prefix}_fold{fold+1}.pt", map_location=device))
        model.eval()

        with torch.no_grad():
            for gene_bodies, promoters, gene_ids in train_loader:
                gene_bodies_cuda = gene_bodies.to(device)
                preds = model(gene_bodies_cuda)
                preds_cpu = preds.cpu()
                for i, gene_id in enumerate(gene_ids):
                    if gene_id not in gene_pred_dict:
                        gene_pred_dict[gene_id] = preds_cpu[i].numpy()
                        gene_body_dict[gene_id] = gene_bodies[i].cpu().numpy()

    avg_promoter = np.stack(list(gene_pred_dict.values()), axis=0).mean(axis=0)
    gene_body_signals = list(gene_body_dict.values())
    max_len = max(arr.shape[0] for arr in gene_body_signals)
    padded_gene_bodies = np.stack([np.pad(arr, (0, max_len - len(arr)), 'constant') for arr in gene_body_signals], axis=0)
    avg_gene_body = padded_gene_bodies.mean(axis=0)

    return avg_promoter, avg_gene_body

# main Execution
if __name__ == "__main__":
    device = torch.device("cuda:2" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Predict for both age groups
    promoter_3yr, gene_body_3yr = predict_and_average(
        "morethan50_nopromoter_nogenebody_cond_3yr_h3k4me3_heartLV.csv",
        "kfold300epochs",
        device
    )
    promoter_34yr, gene_body_34yr = predict_and_average(
        "morethan50_nopromoter_nogenebody_cond_34yr_h3k4me3_heartLV.csv",
        "kfold300epochs_34yr",
        device
    )

    # separate promoter and gene body plots in one frame
    promoter_indices = [1, 2, 3, 4]
    promoter_x = np.arange(2, 6)  # 2, 3, 4, 5
    max_gb = 200
    gene_body_x = np.arange(0, max_gb)

    promoter_y_3yr = promoter_3yr[promoter_indices]
    promoter_y_34yr = promoter_34yr[promoter_indices]
    gene_body_y_3yr = gene_body_3yr[:max_gb]
    gene_body_y_34yr = gene_body_34yr[:max_gb]

    fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharey=False)
    
    plt.rcParams.update({'font.size': 14})

    # Promoter plot
    axs[0].plot(promoter_x, promoter_y_3yr, label="3yr", marker='o', markersize=8, linewidth=3)
    axs[0].plot(promoter_x, promoter_y_34yr, label="34yr", marker='o', markersize=8, linewidth=3)
    axs[0].set_title("Promoter Region (H3K4me3)", fontsize=18, fontweight='bold')
    axs[0].set_xlabel("Position relative to TSS (bins)", fontsize=16)
    axs[0].set_ylabel("Promoter Normalized Signal", fontsize=16)
    axs[0].set_xticks([2, 3, 4, 5])
    axs[0].set_xticklabels(['2', '3', '4', '5'], fontsize=14)
    y_min_promoter = min(np.min(promoter_y_3yr), np.min(promoter_y_34yr)) - 0.01
    y_max_promoter = max(np.max(promoter_y_3yr), np.max(promoter_y_34yr)) + 0.01
    axs[0].set_ylim(y_min_promoter, y_max_promoter)
    axs[0].legend(fontsize=14, frameon=True, fancybox=True, shadow=True)
    axs[0].grid(True, alpha=0.3)
    axs[0].tick_params(axis='both', labelsize=14)

    # Gene body plot
    axs[1].plot(gene_body_x, gene_body_y_3yr, label="3yr", linewidth=2.5)
    axs[1].plot(gene_body_x, gene_body_y_34yr, label="34yr", linewidth=2.5)
    axs[1].set_title("Gene Body Region (H3K4me3)", fontsize=18, fontweight='bold')
    axs[1].set_xlabel("Position relative to TSS (bins)", fontsize=16)
    axs[1].set_ylabel("Gene Body Normalized Signal", fontsize=16)
    x_ticks = np.arange(0, max_gb + 1, 25)
    axs[1].set_xticks(x_ticks)
    axs[1].set_xticklabels([str(int(x)) for x in x_ticks], fontsize=14)
    y_min_genebody = min(np.min(gene_body_y_3yr), np.min(gene_body_y_34yr)) - 0.01
    y_max_genebody = max(np.max(gene_body_y_3yr), np.max(gene_body_y_34yr)) + 0.01
    axs[1].set_ylim(y_min_genebody, y_max_genebody)
    axs[1].legend(fontsize=14, frameon=True, fancybox=True, shadow=True)
    axs[1].grid(True, alpha=0.3)
    axs[1].tick_params(axis='both', labelsize=14)

    plt.tight_layout(pad=2.0)
    plt.savefig("separate_promoter_genebody_signal_h3k4me3_3and34_largefont.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("Plot saved: separate_promoter_genebody_signal_h3k4me3_3and34_largefont.png")

