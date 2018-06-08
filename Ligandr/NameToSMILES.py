"""
NameToSMILES.py
Converts a PDB alphanumeric code to a Canonical Isomeric SMILES 
string with explicit hydrogens at a certain pH

"""

import pypdb
import openbabel


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

def addH(smile_string, pH=7.0):
    """
    Protonates a SMILES string for a given pH

    Parameters
    ------------
    smile_string : str, required
        The SMILES format string of a chemical
    
    pH : double, optional default = 7.0
        The pH at which protonation should occur
    
    Returns
    --------
    protonated_smiles : str
        The SMILES string representing the protonated ligand
    """
    mol = openbabel.OBMol()
    converter = openbabel.OBConversion()
    converter.SetInAndOutFormats("CAN", "CAN")      #Necessary for writing and reading strings
    converter.ReadString(mol,smile_string)
    mol.CorrectForPH(pH)
    mol.AddHydrogens()
    converter.AddOption("h")      #Explicitly add hydrogens
    protonated_smiles = converter.WriteString(mol) 
    return protonated_smiles
