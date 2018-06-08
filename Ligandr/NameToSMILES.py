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
    mol = openbabel.OBMol()
    converter = openbabel.OBConversion()
    converter.SetInAndOutFormats("CAN", "CAN")
    converter.ReadString(mol,smile_string)
    mol.CorrectForPH(pH)
    mol.AddHydrogens()
    print(mol.NumAtoms())
    return ""

def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote

