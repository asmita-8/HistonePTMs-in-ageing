import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
import pandas as pd
import numpy as np
from torch.nn.utils.rnn import pad_sequence
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt
import os

# Focal loss implementation
class FocalLoss(nn.Module):
    def __init__(self, alpha=1, gamma=2, reduction='mean'):
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction

    def forward(self, inputs, targets):
        BCE_loss = F.binary_cross_entropy(inputs, targets, reduction='none')
        pt = torch.where(targets == 1, inputs, 1 - inputs)
        focal_loss = self.alpha * (1 - pt) ** self.gamma * BCE_loss
        return focal_loss.mean() if self.reduction == 'mean' else focal_loss.sum()

# Dataset class
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

# Collate function
def collate_fn(batch):
    sequences, labels = zip(*batch)
    padded_sequences = pad_sequence(sequences, batch_first=True, padding_value=0)
    labels = torch.stack(labels)
    return padded_sequences, labels

# CNN model
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

# Training and testing functions
def train(dataloader, model, loss_function, optimizer):
    model.train()
    total_loss = 0
    correct = 0
    total = 0
    for X, y in dataloader:
        X, y = X.to(device), y.to(device)
        pred = model(X)
        loss = loss_function(pred, y)
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        total_loss += loss.item()
        correct += (pred >= 0.5).float().eq(y).sum().item()
        total += y.numel()
    accuracy = correct / total
    return total_loss / len(dataloader), accuracy

def test(dataloader, model, loss_function):
    model.eval()
    total_loss = 0
    correct = 0
    total = 0
    with torch.no_grad():
        for X, y in dataloader:
            X, y = X.to(device), y.to(device)
            pred = model(X)
            loss = loss_function(pred, y)
            total_loss += loss.item()
            correct += (pred >= 0.5).float().eq(y).sum().item()
            total += y.numel()
    accuracy = correct / total
    return total_loss / len(dataloader), accuracy

# Main K-Fold Training Logic
filepath = "morethan50_nopromoter_nogenebody_cond_34yr_h3k4me3_heartLV.csv"
dataset = GeneDataset(filepath, sep="\t")
k_folds = 5
batch_size = 64
num_epochs = 300
loss_function = FocalLoss()
device = torch.device("cuda:2" if torch.cuda.device_count() > 2 else "cuda" if torch.cuda.is_available() else "cpu")

kf = KFold(n_splits=k_folds, shuffle=True, random_state=42)
all_gene_ids = dataset.gene_ids
fold_accuracies = []

for fold, (train_idx, val_idx) in enumerate(kf.split(all_gene_ids)):
    print(f"\n--- Fold {fold + 1}/{k_folds} ---")
    train_ids = [all_gene_ids[i] for i in train_idx]
    val_ids = [all_gene_ids[i] for i in val_idx]

    train_dataset = GeneDataset(filepath, gene_subset=train_ids)
    val_dataset = GeneDataset(filepath, gene_subset=val_ids)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, collate_fn=collate_fn)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, collate_fn=collate_fn)

    model = H3K4me3CNN().to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=1e-3)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=10, verbose=True)

    best_val_accuracy = 0
    for epoch in range(num_epochs):
        train_loss, train_acc = train(train_loader, model, loss_function, optimizer)
        val_loss, val_acc = test(val_loader, model, loss_function)
        scheduler.step(val_acc)

        if (epoch + 1) % 50 == 0:
            print(f"Epoch {epoch+1}:")
            print(f"  Train Acc = {train_acc:.4f}, Train Loss = {train_loss:.4f}")
            print(f"  Val   Acc = {val_acc:.4f}, Val   Loss = {val_loss:.4f}")

        if val_acc > best_val_accuracy:
            best_val_accuracy = val_acc

    # Save model
    model_save_path = f"kfold300epochs_34yr_fold{fold + 1}.pt"
    torch.save(model.state_dict(), model_save_path)
    print(f"Best Accuracy for Fold {fold + 1}: {best_val_accuracy:.4f}")
    print(f"Model saved to {model_save_path}")

    fold_accuracies.append(best_val_accuracy)

print(f"\nAverage Accuracy across {k_folds} folds: {np.mean(fold_accuracies):.4f}")

