from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from app.services.standardize import standardize_smiles
from app.services.fingerprint_generation import generate_morgan_fingerprint
from app.models.schemas import FingerprintResponse


class StringInput(BaseModel):
    data: str

router = APIRouter()

@router.get("/{SMILES}")
async def get_fingerprint(SMILES:str):
    """
    Generate a molecular fingerprint from a SMILES string.

    Args:
        SMILES string input

    Response:
        Standardized SMILES and Fingerprint as a binary vector
    """
    
    
    print(f"Received data: {SMILES}")
    
    try:
        #step 1: get standardized SMILES
        standardized = standardize_smiles(smiles=SMILES)
        if standardized is None:
            raise ValueError("Invalid SMILES or standardization failed.")

        #step 2: generate fingerprint from standardized SMILES
        fingerprint = generate_morgan_fingerprint(
                smiles=standardized,
                radius=2,
                n_bits=1024
            )
        if fingerprint is None:
            raise ValueError("Fingerprint failed.")

        return FingerprintResponse(standardized_smiles=standardized, fingerprint=fingerprint)

    except Exception as e:
        print(f"Error during fingerprint generation: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Error generating fingerprint: {str(e)}"
        )
