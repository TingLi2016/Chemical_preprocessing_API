from pydantic import BaseModel
from typing import List

class FingerprintResponse(BaseModel):
    """response schema for generating fingerprint."""
    standardized_smiles: str
    fingerprint: List[int]
