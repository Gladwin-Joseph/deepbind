import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool

class DrugGNN(nn.Module):
    """Graph Neural Network for processing drug molecules"""
    def __init__(self, num_node_features, hidden_dim=128):
        super(DrugGNN, self).__init__()
        self.conv1 = GCNConv(num_node_features, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.conv3 = GCNConv(hidden_dim, hidden_dim)
        self.lin = nn.Linear(hidden_dim, hidden_dim)
    
    def forward(self, x, edge_index, batch):
        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = F.relu(self.conv3(x, edge_index))
        x = global_mean_pool(x, batch)
        x = self.lin(x)
        return x

class DTIModel(nn.Module):
    """Drug-Target Interaction prediction model"""
    def __init__(self, drug_in_dim, prot_emb_dim, hidden_dim=128):
        super(DTIModel, self).__init__()
        self.drug_gnn = DrugGNN(drug_in_dim, hidden_dim)
        self.fc1 = nn.Linear(hidden_dim + prot_emb_dim, 256)
        self.fc2 = nn.Linear(256, 128)
        self.out = nn.Linear(128, 1)
    
    def forward(self, data, prot_emb):
        drug_emb = self.drug_gnn(data.x, data.edge_index, data.batch)
        combined = torch.cat((drug_emb, prot_emb), dim=1)
        x = F.relu(self.fc1(combined))
        x = F.dropout(x, p=0.3, training=self.training)
        x = F.relu(self.fc2(x))
        out = self.out(x)
        return out.squeeze()

class Generator(nn.Module):
    """GAN Generator for creating drug embeddings"""
    def __init__(self, latent_dim, output_dim):
        super(Generator, self).__init__()
        self.fc1 = nn.Linear(latent_dim, 256)
        self.fc2 = nn.Linear(256, 512)
        self.fc3 = nn.Linear(512, output_dim)
    
    def forward(self, z):
        x = F.relu(self.fc1(z))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x

class Discriminator(nn.Module):
    """GAN Discriminator"""
    def __init__(self, input_dim):
        super(Discriminator, self).__init__()
        self.fc1 = nn.Linear(input_dim, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 1)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = torch.sigmoid(self.fc3(x))
        return x