from pydantic import BaseModel, Field
from typing import Optional, List

class DTIPredictionRequest(BaseModel):
    """Request schema for DTI prediction"""
    smiles: str = Field(..., description="Drug SMILES string")
    protein_sequence: str = Field(..., description="Protein FASTA sequence")
    
    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "COC1=C(C=C2C(=C1)CCN=C2C3=CC(=C(C=C3)Cl)Cl)Cl",
                "protein_sequence": "MTVKTEAAKGTLTYSRMRGMVAILIAFMKQRRMGLNDFIQKIANNS"
            }
        }

class DTIPredictionResponse(BaseModel):
    """Response schema for DTI prediction"""
    binding_affinity: float = Field(..., description="Predicted binding affinity score")
    smiles: str
    protein_length: int
    
class DrugDiscoveryRequest(BaseModel):
    """Request schema for drug discovery"""
    protein_sequence: str = Field(..., description="Target protein FASTA sequence")
    num_candidates: int = Field(default=10, ge=1, le=100, description="Number of drug candidates to generate")
    
    class Config:
        json_schema_extra = {
            "example": {
                "protein_sequence": "MTVKTEAAKGTLTYSRMRGMVAILIAFMKQRRMGLNDFIQKIANNS",
                "num_candidates": 10
            }
        }

class DrugCandidate(BaseModel):
    """Single drug candidate"""
    rank: int
    smiles: str
    affinity_score: float
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None

class DrugDiscoveryResponse(BaseModel):
    """Response schema for drug discovery"""
    protein_length: int
    total_generated: int
    valid_molecules: int
    candidates: List[DrugCandidate]
    
class HealthResponse(BaseModel):
    """Health check response"""
    status: str
    cuda_available: bool
    device: str