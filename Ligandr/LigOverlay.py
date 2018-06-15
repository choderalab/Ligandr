"""
Places a ligand or multiple ligands in place of another ligand on a protein
"""
import sys
import tempfile
import traceback
from tempfile import _TemporaryFileWrapper
import NameToSMILES

import mdtraj
import numpy as np
from openeye import oechem, oedocking
import pypdb


def dock_molecule(trajectory_frame, original_ligands, new_ligands=None, keep_ions=None) -> mdtraj.Trajectory:
    """
    Docks ligands or a list of ligands in place of another ligand

    Parameters
    ----------
    trajectory_frame: mdtraj.Trajectory, required
        The trajectory to dock to- water will be removed, but every other atoms will stay the same
    """
    # Produces pdb file in memory
    pdb_dummy = tempfile.NamedTemporaryFile(suffix=".pdb")  # type: _TemporaryFileWrapper
    overall_dry_frame = trajectory_frame.remove_solvent(exclude=keep_ions)
    overall_dry_frame.save(pdb_dummy.name)

    # Read pdb to openeye
    complex_prot = oechem.OEGraphMol()
    ins = oechem.oemolistream()
    ins.open(pdb_dummy.name)
    oechem.OEReadMolecule(ins, complex_prot)
    pdb_dummy.close()
    print(complex_prot.NumAtoms())

    # Being processing a ligand
    for ligand in original_ligands:
        try:
            ligand_traj = overall_dry_frame.atom_slice(overall_dry_frame.top.select(ligand))  # type: mdtraj.Trajectory
        except IndexError as e:
            print("Ligand selection was invalid")
            print("Exception: %s" % str(e), file=sys.stderr)
            traceback.print_tb(e.__traceback__, file=sys.stderr)
            sys.exit(1)
        if ligand_traj.n_residues > 1:
            raise IndexError("More than one ligand residue found. Only one expected")
        # Find approximate location of binding pocket
        print(len(ligand_traj.xyz[0]))
        ligand_center = np.average(ligand_traj.xyz[0], axis=0) * 10

        # Explicit value necessary for generated methods
        x = float(ligand_center[0])
        y = float(ligand_center[1])
        z = float(ligand_center[2])

        # Create receptor molecule
        receptor = oechem.OEGraphMol()
        oedocking.OEMakeReceptor(receptor, complex_prot, x, y, z)

        # Write out and open ligand in OpenEye
        ligand_dummy = tempfile.NamedTemporaryFile(suffix=".pdb")
        ligand_traj.save(ligand_dummy.name)
        lig_initial = oechem.OEGraphMol()
        ins.open(ligand_dummy.name)
        oechem.OEReadMolecule(ins, lig_initial)
        oechem.OEAddExplicitHydrogens(lig_initial)

        # Validate if ligand is correct by checking against SMILES
        lig_SMILES = NameToSMILES.getSMILE(ligand_traj.top.residue(0).name)
        neutral_ph_smiles = NameToSMILES.addH(lig_SMILES)
        validation_ligand = oechem.OEGraphMol()
        oechem.OESmilesToMol(validation_ligand, neutral_ph_smiles)
        if oechem.OEMolToSmiles(validation_ligand) != oechem.OEMolToSmiles(lig_initial):
            __docking_internal(receptor, lig_initial, validation_ligand)


def __docking_internal(receptor, bound_ligand, docking_ligand):
    print("End")
