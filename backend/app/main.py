from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import torch
import os
from pathlib import Path

from .schemas import (
    DTIPredictionRequest,
    DTIPredictionResponse,
    DrugDiscoveryRequest,
    DrugDiscoveryResponse,
    DrugCandidate,
    HealthResponse
)
from .models import DTIModel, Generator
from .config import DEVICE, LOAD_MODELS_ON_STARTUP, DTI_MODEL_PATH, GENERATOR_MODEL_PATH

# Initialize FastAPI app
app = FastAPI(
    title="Drug Discovery API",
    description="AI-powered drug discovery platform",
    version="1.0.0"
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # For Render deployment
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global variables - Lazy loaded
dti_predictor = None
drug_generator = None
models_loaded = False

def load_models():
    """Lazy load models only when needed"""
    global dti_predictor, drug_generator, models_loaded
    
    if models_loaded:
        return True
    
    try:
        print("🔄 Loading models...")
        
        # Check if model files exist
        if not os.path.exists(DTI_MODEL_PATH):
            print(f"⚠️ DTI model not found at {DTI_MODEL_PATH}")
            return False
        
        if not os.path.exists(GENERATOR_MODEL_PATH):
            print(f"⚠️ Generator model not found at {GENERATOR_MODEL_PATH}")
            return False
        
        # Load DTI model
        from .utils.dti_predictor import DTIPredictor
        dti_model = DTIModel(drug_in_dim=1, prot_emb_dim=1024)
        dti_model.load_state_dict(torch.load(DTI_MODEL_PATH, map_location=DEVICE))
        dti_predictor = DTIPredictor(dti_model, DEVICE)
        print("✅ DTI model loaded")
        
        # Load Generator
        from .utils.drug_generator import DrugGenerator
        generator = Generator(latent_dim=56, output_dim=128)
        generator.load_state_dict(torch.load(GENERATOR_MODEL_PATH, map_location=DEVICE))
        drug_generator = DrugGenerator(generator, dti_model, DEVICE)
        print("✅ Generator model loaded")
        
        models_loaded = True
        return True
        
    except Exception as e:
        print(f"❌ Error loading models: {str(e)}")
        return False

@app.on_event("startup")
async def startup_event():
    """Load models on startup if flag is set"""
    print("=" * 60)
    print("INITIALIZING DRUG DISCOVERY API")
    print("=" * 60)
    print(f"✓ Device: {DEVICE}")
    
    if LOAD_MODELS_ON_STARTUP:
        load_models()
    else:
        print("⚠️ Models will be loaded on first request (lazy loading)")
    
    print("=" * 60)
    print("API READY!")
    print("=" * 60)

@app.get("/", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    return HealthResponse(
        status="healthy",
        cuda_available=torch.cuda.is_available(),
        device=str(DEVICE)
    )

@app.post("/api/predict-dti", response_model=DTIPredictionResponse)
async def predict_dti(request: DTIPredictionRequest):
    """Predict drug-target interaction binding affinity"""
    
    # Lazy load models if not already loaded
    if not models_loaded:
        success = load_models()
        if not success:
            raise HTTPException(
                status_code=503, 
                detail="Models not available. Please contact administrator."
            )
    
    try:
        print(f"\n→ DTI Prediction Request")
        
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
    """Generate novel drug candidates for a target protein"""
    
    # Lazy load models if not already loaded
    if not models_loaded:
        success = load_models()
        if not success:
            raise HTTPException(
                status_code=503, 
                detail="Models not available. Please contact administrator."
            )
    
    try:
        print(f"\n→ Drug Discovery Request")
        
        candidates = drug_generator.generate_candidates(
            request.protein_sequence,
            request.num_candidates
        )
        
        print(f"✓ Generated {len(candidates)} valid candidates")
        
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
        "device": str(DEVICE),
        "models_loaded": models_loaded,
        "dti_model_loaded": dti_predictor is not None,
        "generator_loaded": drug_generator is not None
    }

# Include render router if exists
try:
    from .render import router as render_router
    app.include_router(render_router)
    print("✓ Render router included")
except ImportError:
    print("⚠️ Render router not available")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)