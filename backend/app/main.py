from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import torch
import os
from pathlib import Path
from .render import router as render_router

from .schemas import (
    DTIPredictionRequest,
    DTIPredictionResponse,
    DrugDiscoveryRequest,
    DrugDiscoveryResponse,
    DrugCandidate,
    HealthResponse
)
from .models import DTIModel, Generator, Discriminator
from .utils.dti_predictor import DTIPredictor
from .utils.drug_generator import DrugGenerator

# Initialize FastAPI app
app = FastAPI(
    title="Drug Discovery API",
    description="AI-powered drug discovery platform with DTI prediction and drug generation",
    version="1.0.0"
)

app.include_router(render_router)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000","http://localhost:5173","https://drug-discovery-frontend.onrender.com"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global variables
device = None
dti_predictor = None
drug_generator = None

@app.on_event("startup")
async def startup_event():
    """Load models on startup"""
    global device, dti_predictor, drug_generator
    
    print("=" * 60)
    print("INITIALIZING DRUG DISCOVERY API")
    print("=" * 60)
    
    # Set device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"✓ Device: {device}")
    
    if torch.cuda.is_available():
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
        print(f"  CUDA Version: {torch.version.cuda}")
    
    # Model paths
    base_path = Path(__file__).parent.parent / "models"
    dti_model_path = base_path / "best_dti_model.pth"
    generator_path = base_path / "generator_model.pth"
    
    # Check if models exist
    if not dti_model_path.exists():
        raise FileNotFoundError(f"DTI model not found at {dti_model_path}")
    if not generator_path.exists():
        raise FileNotFoundError(f"Generator model not found at {generator_path}")
    
    print(f"\n→ Loading DTI model from: {dti_model_path}")
    
    # Load DTI model
    dti_model = DTIModel(drug_in_dim=1, prot_emb_dim=1024)
    dti_model.load_state_dict(torch.load(dti_model_path, map_location=device))
    dti_predictor = DTIPredictor(dti_model, device)
    print("✓ DTI model loaded successfully")
    
    print(f"\n→ Loading Generator model from: {generator_path}")
    
    # Load Generator
    generator = Generator(latent_dim=56, output_dim=128)
    generator.load_state_dict(torch.load(generator_path, map_location=device))
    drug_generator = DrugGenerator(generator, dti_model, device)
    print("✓ Generator model loaded successfully")
    
    print("\n" + "=" * 60)
    print("API READY!")
    print("=" * 60)

@app.get("/", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    return HealthResponse(
        status="healthy",
        cuda_available=torch.cuda.is_available(),
        device=str(device)
    )

@app.post("/api/predict-dti", response_model=DTIPredictionResponse)
async def predict_dti(request: DTIPredictionRequest):
    """
    Predict drug-target interaction binding affinity
    
    - **smiles**: Drug molecule in SMILES format
    - **protein_sequence**: Target protein amino acid sequence
    """
    try:
        print(f"\n→ DTI Prediction Request")
        print(f"  SMILES: {request.smiles}")
        print(f"  Protein length: {len(request.protein_sequence)}")
        
        affinity = dti_predictor.predict(
            request.smiles,
            request.protein_sequence
        )
        
        print(f"✓ Predicted affinity: {affinity:.4f}")
        
        return DTIPredictionResponse(
            binding_affinity=affinity,
            smiles=request.smiles,
            protein_length=len(request.protein_sequence)
        )
        
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction failed: {str(e)}")

@app.post("/api/discover-drugs", response_model=DrugDiscoveryResponse)
async def discover_drugs(request: DrugDiscoveryRequest):
    """
    Generate novel drug candidates for a target protein
    
    - **protein_sequence**: Target protein amino acid sequence
    - **num_candidates**: Number of drug candidates to generate (1-100)
    """
    try:
        print(f"\n→ Drug Discovery Request")
        print(f"  Protein length: {len(request.protein_sequence)}")
        print(f"  Candidates requested: {request.num_candidates}")
        
        candidates = drug_generator.generate_candidates(
            request.protein_sequence,
            request.num_candidates
        )
        
        print(f"✓ Generated {len(candidates)} valid candidates")
        
        # Format response
        drug_candidates = [
            DrugCandidate(
                rank=idx + 1,
                smiles=c['smiles'],
                affinity_score=c['affinity'],
                molecular_weight=c['mw'],
                logp=c['logp']
            )
            for idx, c in enumerate(candidates)
        ]
        
        return DrugDiscoveryResponse(
            protein_length=len(request.protein_sequence),
            total_generated=request.num_candidates * 3,
            valid_molecules=len(candidates),
            candidates=drug_candidates
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Drug generation failed: {str(e)}")

@app.get("/api/health")
async def detailed_health():
    """Detailed health check with model status"""
    return {
        "status": "healthy",
        "cuda_available": torch.cuda.is_available(),
        "device": str(device),
        "gpu_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "dti_model_loaded": dti_predictor is not None,
        "generator_loaded": drug_generator is not None
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)