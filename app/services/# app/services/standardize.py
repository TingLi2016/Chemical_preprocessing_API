# app/services/standardize.py

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import Optional

# Suppress RDKit warnings and errors globally
RDLogger.DisableLog('rdApp.*')


def standardize_smiles(smiles: str) -> Optional[str]:
    """
    Standardizes a chemical structure from a SMILES string.

    Steps include:
    - Sanitization and kekulization check
    - Fragment removal
    - Functional group normalization
    - Reionization
    - Tautomer canonicalization

    Args:
        smiles (str): Input SMILES string.

    Returns:
        Optional[str]: Standardized SMILES string or None if failed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[Standardization] Invalid SMILES input: {smiles}")
        return None

    # Pre-processing: Sanitize and kekulize
    try:
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        print(f"[Standardization] Pre-processing error: {e}")
        return None

    try:
        # Step 1: Remove fragments (keep largest)
        fragment_remover = rdMolStandardize.FragmentRemover()
        mol = fragment_remover.remove(mol)

        # Step 2: Normalize functional groups
        normalizer = rdMolStandardize.Normalizer()
        mol = normalizer.normalize(mol)

        # Step 3: Reionization
        reionizer = rdMolStandardize.Reionizer()
        mol = reionizer.reionize(mol)

        # Step 4: Tautomer canonicalization
        tautomer_enum = rdMolStandardize.TautomerEnumerator()
        mol = tautomer_enum.Canonicalize(mol)

        # Output standardized SMILES
        return Chem.MolToSmiles(mol)

    except Exception as e:
        print(f"[Standardization] Error during standardization: {e}")
        return None
