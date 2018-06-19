"""
Places a ligand or multiple ligands in place of another ligand on a protein
"""
import sys
import tempfile
import traceback
from typing import List

import NameToSMILES

import mdtraj
import numpy as np
from mdtraj import Trajectory
from openeye import oechem, oedocking, oeomega, oeff
from openeye.oechem import OEMol, OEGraphMol
import itertools
import subprocess


def dock_molecule(trajectory_frame, original_ligands, new_ligands, title_list, keep_ions=None,
                  override_replacement=False) -> List[
    mdtraj.Trajectory]:
    """
    Docks ligands or a list of ligands in place of another ligand

    Parameters
    ----------
    trajectory_frame: mdtraj.Trajectory, required
        The trajectory to dock to- water will be removed, but every other atoms will stay the same

    original_ligands: [str], required
        The original ligand to replace, must be an mdtraj selection string. If more than
        one residue or no residue is selected, program will throw an index error

    new_ligands: [[str]], required
        For each original ligand, a list of SMILES string representing the new ligand to dock.

    title_list: [[str]], required
        For each new ligand, the title the user wishes to have.

    keep_ions: [str], optional, default:None
        A list of ions to keep when removing solvent. Must be an mdtraj selection.

    override_replacement: bool, optional, default:False
        Setting to true will prevent replacing ligand with canonical molecule from PDB.

    Returns
    -------
    output: [mdtraj.Trajectory]
        A list of trajectories for all combinations of ligands
    """
    # Check inputs
    if len(original_ligands) != len(new_ligands):
        raise ValueError("Original and new ligand arrays must have the same length.")
    # Produces pdb file in memory
    overall_dry_frame = trajectory_frame.remove_solvent(exclude=keep_ions)  # type: mdtraj.Trajectory
    complex_prot = __trajectory_to_oemol__(overall_dry_frame)

    # Bookkeeping
    ligand_stream = []
    output = []
    i = 0
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
        ligand_center = np.average(ligand_traj.xyz[0], axis=0) * 10

        # Explicit value necessary for generated methods
        x = float(ligand_center[0])
        y = float(ligand_center[1])
        z = float(ligand_center[2])

        # Create receptor molecule
        receptor = oechem.OEGraphMol()
        oedocking.OEMakeReceptor(receptor, complex_prot, x, y, z)

        # Write out and open ligand in OpenEye
        lig_initial = __trajectory_to_oemol__(ligand_traj)
        lig_initial.SetTitle(ligand_traj.top.residue(0).name)


        # Validate if ligand is correct by checking against SMILES
        lig_smiles = NameToSMILES.getSMILE(ligand_traj.top.residue(0).name)
        if lig_smiles is not None:
            neutral_ph_smiles = NameToSMILES.addH(lig_smiles)
            validation_ligand = oechem.OEGraphMol()
            validation_ligand.SetTitle(ligand_traj.top.residue(0).name)
            oechem.OESmilesToMol(validation_ligand, neutral_ph_smiles)
            if (oechem.OEMolToSmiles(validation_ligand) != oechem.OEMolToSmiles(
                    lig_initial)) and not override_replacement:
                ligand_posed = __docking_internal(receptor, lig_initial, validation_ligand)
                complex_prot = __trajectory_to_oemol__(
                    overall_dry_frame.atom_slice(overall_dry_frame.top.select("not ( %s )" % ligand)))
                oechem.OEAddMols(receptor, ligand_posed)
                print(
                    "Warning: Original ligand had either incorrect stereochemistry or missing heavy atoms. Ligand from pdb database placed instead. ")

            else:
                ligand_posed = lig_initial
                print("stable")
        else:
            print("ligand not found in PDB: No validation and continuing")
        # Dock new ligands
        ligand_position = []  # type: List[OEGraphMol]
        j = 0
        for new_ligand in new_ligands[i]:
            new_structure = oechem.OEGraphMol()
            oechem.OESmilesToMol(new_structure, new_ligand)
            new_structure.SetTitle("ligand replacement")
            docked_ligand = __docking_internal(receptor, ligand_posed, new_structure)
            docked_ligand.SetTitle(title_list[i][j])
            ligand_position.append(docked_ligand)
            j = j + 1
        ligand_stream.append(ligand_position)
        overall_dry_frame.atom_slice(overall_dry_frame.top.select("not ( %s )" % ligand), inplace=True)
        i = i + 1
    # Produce output trajectories for all combos of ligands
    for ligand_combo in itertools.product(*ligand_stream):
        new_frame_prot = overall_dry_frame  # type: Trajectory
        for built_ligand in ligand_combo:
            # Rename ligand
            ligand_dummy = tempfile.NamedTemporaryFile(suffix=".pdb")
            outputter_lig = oechem.oemolostream(ligand_dummy.name)
            oechem.OEWriteMolecule(outputter_lig, built_ligand)
            ligand_final_traj = mdtraj.load(ligand_dummy.name)
            ligand_dummy.close()
            for residue in ligand_final_traj.top.residues:
                residue.name = built_ligand.GetTitle()
            new_frame_prot = new_frame_prot.stack(ligand_final_traj)
        output.append(new_frame_prot)
    print("End")
    return output


def __docking_internal(receptor: oechem.OEGraphMol, bound_ligand: oechem.OEGraphMol, docking_ligand: oechem.OEGraphMol):
    """
    Internal method for docking to a receptor given a bound ligand and ligand to dock to

    Parameters
    ----------
    receptor: oechem.OEGraphMol, required
        An oechem receptor to dock to

    bound_ligand: oechem.OEGraphMol, required
        The ligand bound to the current receptor pocket

    docking_ligand: oechem.OEGraphMol, required
        The ligand to dock.

    Returns
    -------
    A oechem molecules with coordinates of its docking

    """
    # Create posit config
    oedocking.OEReceptorSetBoundLigand(receptor, bound_ligand)
    poser = oedocking.OEPosit()
    poser.Initialize(receptor)

    # Create multiple conformations
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfs(1000)
    omega_driver = oeomega.OEOmega(omegaOpts)
    conformer_docking = oechem.OEMol()  # type: OEMol
    if not omega_driver(conformer_docking):
        single_sdf = tempfile.NamedTemporaryFile(suffix=".sdf")
        double_sdf = tempfile.NamedTemporaryFile(suffix=".sdf")
        smiles = NameToSMILES.addH(oechem.OEMolToSmiles(docking_ligand))
        subprocess.run("obabel -:'%s' -O %s --gen3d; obabel %s -O %s --confab --conf = 100" % (
            smiles, single_sdf.name, single_sdf.name, double_sdf.name), shell=True)
        ins = oechem.oemolistream()
        ins.SetConfTest(oechem.OEOmegaConfTest())
        ins.open(double_sdf.name)
        oechem.OEReadMolecule(ins, conformer_docking)

    # Failed conformer generation
    if conformer_docking.NumConfs() <= 1:
        raise TypeError('Conformer generation failed. Only %s conformer(s) present' % conformer_docking.NumConfs())
    # Dock and get top conformer
    posed_ligand = oechem.OEGraphMol()  # type: OEGraphMol
    poser.DockMultiConformerMolecule(posed_ligand, conformer_docking)
    return posed_ligand


def __trajectory_to_oemol__(trajectory):
    """
    Converts an mdtrajectory into an openeye molecule

    Parameters
    ---------
    trajectory: mdtraj.Trajectory, required
        The trajectory to return

    Returns
    -------
    mol: oechem.OEGraphMol
        The molecule that is returned
    """
    dummy_file = tempfile.NamedTemporaryFile(suffix=".pdb")
    trajectory.save(dummy_file.name)
    mol = oechem.OEGraphMol()
    ins = oechem.oemolistream()
    ins.open(dummy_file.name)
    oechem.OEReadMolecule(ins, mol)
    if not mol.NumAtoms() >= 1:
        raise AssertionError("Error reading in molecule")
    dummy_file.close()
    return mol
