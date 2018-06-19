"""
NameToSMILES.py
Converts a PDB alphanumeric code to a Canonical Isomeric SMILES 
string with explicit hydrogens at pH 7.4

"""

from openeye import oechem, oequacpac
import requests


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
    url_string = "https://www.rcsb.org/ligand/" + chemName
    response = requests.get(url_string)
    long_string = response.content.decode()
    containing_string = long_string[long_string.find("Isomeric SMILES"): long_string.find("InChI")]
    suffixed_string = containing_string[54:]
    smiles_string = suffixed_string[:suffixed_string.find("<")]
    return smiles_string


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
    non_protonated = oechem.OEMolToSmiles(ligand)
    return non_protonated
