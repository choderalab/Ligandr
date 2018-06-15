"""
NameToSMILES.py
Converts a PDB alphanumeric code to a Canonical Isomeric SMILES 
string with explicit hydrogens at pH 7.4

"""

import pypdb
from openeye import oechem, oequacpac
from rdkit import Chem



def getSMILE(chemName):
    """
    Finds the SMILES from a 3 or 4 letter code of a chemical located in the PDB
    
    Parameters
    ----------
    chemName : str, required
        The string representing a chemical in the PDB database

    Returns
    --------
    smiles : str
        SMILES representation of chemical 
    """
    chem_detail = pypdb.describe_chemical(chemName)
    smiles = chem_detail["describeHet"]["ligandInfo"]["ligand"]["smiles"]
    return smiles


def addH(smile_string):
    """
    Protonates a SMILES string for a pH at 7.4

    Parameters
    ------------
    smile_string : str, required
        The SMILES format string of a chemical
    Returns
    --------
    protonated_smiles : str
        The SMILES string representing the protonated ligand
    """
    ligand = oechem.OEGraphMol()
    oechem.OESmilesToMol(ligand, smile_string)
    oequacpac.OESetNeutralpHModel(ligand)
    non_protonated = oechem.OECreateIsoSmiString(ligand)
    ligand_ph = Chem.MolFromSmiles(non_protonated)
    protonated_smiles = Chem.MolToSmiles(ligand_ph, allHsExplicit=True)
    return protonated_smiles
