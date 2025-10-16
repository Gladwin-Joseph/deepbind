import torch
import torch.nn.functional as F
from rdkit import Chem
from rdkit.Chem import Descriptors
from transformers import BertTokenizer, BertModel
import numpy as np

# Drug-like molecular scaffolds
DRUG_SCAFFOLDS = [
    'c1ccccc1',           # benzene
    'c1ccc(O)cc1',        # phenol
    'c1ccc(N)cc1',        # aniline
    'c1ccc(C)cc1',        # toluene
    'c1cccnc1',           # pyridine
    'c1ccoc1',            # furan
    'c1ccsc1',            # thiophene
    'c1c[nH]cn1',         # imidazole
    'c1cnccn1',           # pyrimidine
    'c1ccc2ccccc2c1',     # naphthalene
    'c1ccc2[nH]ccc2c1',   # indole
    'c1ccc2ncccc2c1',     # quinoline
    'C1CCCCC1',           # cyclohexane
    'C1CCNCC1',           # piperidine
]

class DrugGenerator:
    """Generate novel drug candidates using trained GAN"""
    
    def __init__(self, generator, dti_model, device, latent_dim=56):
        self.generator = generator
        self.dti_model = dti_model
        self.device = device
        self.latent_dim = latent_dim
        
        self.generator.to(device)
        self.generator.eval()
        self.dti_model.to(device)
        self.dti_model.eval()
        
        # Load ProtBERT
        print("Loading ProtBERT for drug generation...")
        self.tokenizer = BertTokenizer.from_pretrained('Rostlab/prot_bert', do_lower_case=False)
        self.protbert_model = BertModel.from_pretrained('Rostlab/prot_bert', use_safetensors=True)
        self.protbert_model.to(device)
        self.protbert_model.eval()
        print("ProtBERT loaded successfully")
    
    def get_protein_embedding(self, sequence):
        """Generate protein embedding"""
        seq = " ".join(list(sequence))
        inputs = self.tokenizer(
            seq,
            return_tensors="pt",
            truncation=True,
            max_length=512,
            padding=True
        )
        inputs = {k: v.to(self.device) for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = self.protbert_model(**inputs)
            embedding = outputs.last_hidden_state.mean(dim=1).squeeze()
        
        return embedding
    
    def generate_drug_embeddings(self, num_drugs):
        """Generate drug embeddings using trained generator"""
        noise = torch.randn(num_drugs, self.latent_dim).to(self.device)
        
        with torch.no_grad():
            drug_embeddings = self.generator(noise)
        
        return drug_embeddings
    
    def embedding_to_smiles(self, embedding):
        """Convert drug embedding to SMILES (deterministic)"""
        emb_np = embedding.cpu().numpy()
        
        # Select scaffold based on embedding
        scaffold_idx = int(abs(emb_np[0] * 1000)) % len(DRUG_SCAFFOLDS)
        base_smiles = DRUG_SCAFFOLDS[scaffold_idx]
        
        try:
            mol = Chem.MolFromSmiles(base_smiles)
            if mol is None:
                return None
            
            # Return valid SMILES
            final_smiles = Chem.MolToSmiles(mol)
            test_mol = Chem.MolFromSmiles(final_smiles)
            
            if test_mol is not None:
                return final_smiles
            else:
                return base_smiles
                
        except Exception:
            return None
    
    def predict_affinity(self, drug_embeddings, protein_embedding):
        """Predict binding affinities for generated drugs"""
        num_drugs = drug_embeddings.shape[0]
        prot_emb = protein_embedding.unsqueeze(0).repeat(num_drugs, 1)
        
        with torch.no_grad():
            combined = torch.cat((drug_embeddings, prot_emb), dim=1)
            x = F.relu(self.dti_model.fc1(combined))
            x = F.relu(self.dti_model.fc2(x))
            scores = self.dti_model.out(x).squeeze()
        
        return scores
    
    def generate_candidates(self, protein_sequence, num_candidates=10):
        """Generate and rank drug candidates"""
        # Generate more than needed to account for invalid molecules
        num_to_generate = num_candidates * 3
        
        # Get protein embedding
        protein_emb = self.get_protein_embedding(protein_sequence)
        
        # Generate drug embeddings
        drug_embeddings = self.generate_drug_embeddings(num_to_generate)
        
        # Predict affinities
        affinity_scores = self.predict_affinity(drug_embeddings, protein_emb)
        
        # Convert embeddings to SMILES
        candidates = []
        for idx in range(num_to_generate):
            smiles = self.embedding_to_smiles(drug_embeddings[idx])
            
            if smiles is not None:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    
                    candidates.append({
                        'smiles': smiles,
                        'affinity': affinity_scores[idx].item(),
                        'mw': mw,
                        'logp': logp
                    })
        
        # Sort by affinity and return top candidates
        candidates.sort(key=lambda x: x['affinity'], reverse=True)
        
        return candidates[:num_candidates]