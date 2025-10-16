from fastapi import APIRouter, HTTPException
from fastapi.responses import Response
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64

router = APIRouter()

class MoleculeRequest(BaseModel):
    smiles: str
    width: int = 200
    height: int = 200

@router.post("/api/render-molecule")
async def render_molecule(request: MoleculeRequest):
    """Render 2D molecule structure from SMILES"""
    try:
        mol = Chem.MolFromSmiles(request.smiles)
        
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # Generate 2D coordinates
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
        
        # Draw molecule
        img = Draw.MolToImage(mol, size=(request.width, request.height))
        
        # Convert to PNG bytes
        buf = BytesIO()
        img.save(buf, format='PNG')
        buf.seek(0)
        
        return Response(content=buf.getvalue(), media_type="image/png")
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Rendering failed: {str(e)}")

@router.post("/api/render-molecule-base64")
async def render_molecule_base64(request: MoleculeRequest):
    """Render 2D molecule and return as base64"""
    try:
        mol = Chem.MolFromSmiles(request.smiles)
        
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
        
        img = Draw.MolToImage(mol, size=(request.width, request.height))
        
        buf = BytesIO()
        img.save(buf, format='PNG')
        buf.seek(0)
        
        # Convert to base64
        img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
        
        return {"image": f"data:image/png;base64,{img_base64}"}
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Rendering failed: {str(e)}")