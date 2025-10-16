import torch
import torch.nn.functional as F
from torch_geometric.data import Data as GeoData
from rdkit import Chem
from rdkit.Chem import rdchem
from transformers import BertTokenizer, BertModel
import os

class DTIPredictor:
    """Drug-Target Interaction predictor"""
    
    def __init__(self, model, device):
        self.model = model
        self.device = device
        self.model.to(device)
        self.model.eval()
        
        # Load ProtBERT for protein embeddings
        print("Loading ProtBERT model...")
        self.tokenizer = BertTokenizer.from_pretrained('Rostlab/prot_bert', do_lower_case=False)
        self.protbert_model = BertModel.from_pretrained('Rostlab/prot_bert', use_safetensors=True)
        self.protbert_model.to(device)
        self.protbert_model.eval()
        print("ProtBERT model loaded successfully")
    
    def mol_to_graph(self, smiles):
        """Convert SMILES to PyG graph"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Node features (atomic numbers)
        atom_features_list = [[atom.GetAtomicNum()] for atom in mol.GetAtoms()]
        x = torch.tensor(atom_features_list, dtype=torch.float)
        
        # Edge indices and features
        edge_index, edge_attr = [], []
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            edge_index.append([i, j])
            edge_index.append([j, i])
            
            bt = bond.GetBondType()
            bt_num = {
                rdchem.BondType.SINGLE: 1,
                rdchem.BondType.DOUBLE: 2,
                rdchem.BondType.TRIPLE: 3,
                rdchem.BondType.AROMATIC: 4
            }.get(bt, 0)
            
            edge_attr.append([bt_num])
            edge_attr.append([bt_num])
        
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)
        
        return GeoData(x=x, edge_index=edge_index, edge_attr=edge_attr)
    
    def get_protein_embedding(self, sequence):
        """Generate protein embedding using ProtBERT"""
        # Add spaces between amino acids
        seq = " ".join(list(sequence))
        
        # Tokenize and encode
        inputs = self.tokenizer(
            seq, 
            return_tensors="pt", 
            truncation=True, 
            max_length=512,
            padding=True
        )
        inputs = {k: v.to(self.device) for k, v in inputs.items()}
        
        # Get embeddings
        with torch.no_grad():
            outputs = self.protbert_model(**inputs)
            embedding = outputs.last_hidden_state.mean(dim=1).squeeze()
        
        return embedding
    
    def predict(self, smiles, protein_sequence):
        """Predict binding affinity"""
        try:
            # Convert SMILES to graph
            graph = self.mol_to_graph(smiles)
            graph = graph.to(self.device)
            
            # Get protein embedding
            prot_emb = self.get_protein_embedding(protein_sequence)
            prot_emb = prot_emb.unsqueeze(0)  # Add batch dimension
            
            # Create batch tensor for single graph
            batch = torch.zeros(graph.x.shape[0], dtype=torch.long).to(self.device)
            
            # Get drug embedding from GNN
            drug_emb = self.model.drug_gnn(graph.x, graph.edge_index, batch)
            
            # Combine and predict
            combined = torch.cat((drug_emb, prot_emb), dim=1)
            x = F.relu(self.model.fc1(combined))
            x = F.relu(self.model.fc2(x))
            affinity = self.model.out(x).squeeze()
            
            return affinity.item()
            
        except Exception as e:
            raise Exception(f"Prediction failed: {str(e)}")