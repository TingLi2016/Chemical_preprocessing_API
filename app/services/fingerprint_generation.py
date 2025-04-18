# app/services/fingerprint.py

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import numpy as np
from typing import Optional, List

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')


def generate_morgan_fingerprint(smiles: str, radius: int = 2, n_bits: int = 1024) -> Optional[List[int]]:
    """
    Generates a Morgan fingerprint (ECFP) from a SMILES string.

    Args:
        smiles (str): The input SMILES string.
        radius (int, optional): Radius for fingerprint generation. Defaults to 2.
        n_bits (int, optional): Length of the fingerprint bit vector. Defaults to 1024.

    Returns:
        Optional[List[int]]: A list of 0s and 1s representing the fingerprint, or None if SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[Fingerprint] Invalid SMILES: {smiles}")
        return None

    try:
        # Generate Morgan fingerprint as bit vector
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
        arr = np.zeros((n_bits,), dtype=int)
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
        return arr.tolist()

    except Exception as e:
        print(f"[Fingerprint] Error generating fingerprint: {e}")
        return None
