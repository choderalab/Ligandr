import mdtraj
from pdbfixer import PDBFixer
import pdbfixer
from simtk.openmm import app
import tempfile
import ChiralityCheck

ONE_THREE_CODE = \
    {'A': 'ALA',
     'C': 'CYS',
     'D': 'ASP',
     'E': 'GLU',
     'F': 'PHE',
     'G': 'GLY',
     'H': 'HIS',
     'I': 'ILE',
     'K': 'LYS',
     'L': 'LEU',
     'M': 'MET',
     'N': 'ASN',
     'P': 'PRO',
     'Q': 'GLN',
     'R': 'ARG',
     'S': 'SER',
     'T': 'THR',
     'V': 'VAL',
     'W': 'TRP',
     'Y': 'TYR'}


def fix_peptide(pdb_file, seq_dict, pH=7.4, remove_water=True, remove_small_mols=True):
    global ONE_THREE_CODE
    fixer = PDBFixer(filename=pdb_file)
    fixer.sequences.clear()
    for chain in fixer.topology.chains():
        seq = pdbfixer.pdbfixer.Sequence(chain.id, [r.name for r in list(chain.residues())])
        fixer.sequences.append(seq)
    if remove_small_mols:
        fixer.removeHeterogens(not remove_water)
    delete_chains = []
    # Convert single AA codes to three letter code
    for key, value in seq_dict.items():
        if not value or value is None:
            delete_chains.append(key)
        else:
            three_letter = []
            for item in value:
                three_letter.append(ONE_THREE_CODE[item])
            seq_dict[key] = three_letter

    for chain in fixer.topology.chains():
        if chain.index in seq_dict:
            if seq_dict[chain.index] is not None:
                fixer.sequences[chain.index].residues = seq_dict[chain.index]
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)
    fixer.removeChains(delete_chains)
    dummy = tempfile.NamedTemporaryFile(suffix=".pdb")
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open(dummy.name, 'w'))
    product = mdtraj.load(dummy.name)
    problem_cis = ChiralityCheck.check_cispeptide_bond(product)
    problem_chiral = ChiralityCheck.check_chirality(product)
    print("The following problems have been detected:")
    print(problem_cis)
    print(problem_chiral)
    print("Either rerun or find a tool to solve. Perhaps VMD?")
    return product
