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
from openeye import oechem, oedocking, oeomega, oeff
from openeye.oechem import OEMol, OEGraphMol
import itertools
import subprocess


def dock_molecule(trajectory_frame, original_ligands, new_ligands, title_list, keep_ions=None) -> List[
    mdtraj.Trajectory]:
    """
    Docks ligands or a list of ligands in place of another ligand

    Parameters
    ----------
    trajectory_frame: mdtraj.Trajectory, required
        The trajectory to dock to- water will be removed, but every other atoms will stay the same

    original_ligands: str, required
        The original ligand to replace, must be an md traj selection string. If more than
        one residue or no residue is selected, program will throw an index error
    """
    # Produces pdb file in memory
    pdb_dummy = tempfile.NamedTemporaryFile(suffix=".pdb")
    overall_dry_frame = trajectory_frame.remove_solvent(exclude=keep_ions)
    overall_dry_frame.save(pdb_dummy.name)

    # Read pdb to openeye
    complex_prot = oechem.OEGraphMol()
    ins = oechem.oemolistream()
    ins.open(pdb_dummy.name)
    oechem.OEReadMolecule(ins, complex_prot)
    pdb_dummy.close()
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
        ligand_dummy = tempfile.NamedTemporaryFile(suffix=".pdb")
        ligand_traj.save(ligand_dummy.name)
        lig_initial = oechem.OEGraphMol()
        ins.open(ligand_dummy.name)
        oechem.OEReadMolecule(ins, lig_initial)
        oechem.OEAddExplicitHydrogens(lig_initial)
        lig_initial.SetTitle(ligand_traj.top.residue(0).name)
        ligand_dummy.close()

        # Validate if ligand is correct by checking against SMILES
        lig_smiles = NameToSMILES.getSMILE(ligand_traj.top.residue(0).name)
        neutral_ph_smiles = NameToSMILES.addH(lig_smiles)
        validation_ligand = oechem.OEGraphMol()
        validation_ligand.SetTitle(ligand_traj.top.residue(0).name)
        oechem.OESmilesToMol(validation_ligand, neutral_ph_smiles)
        if oechem.OEMolToSmiles(validation_ligand) != oechem.OEMolToSmiles(lig_initial):
            ligand_posed = __docking_internal(receptor, lig_initial, validation_ligand)
            complex_prot = __eraseLigand(receptor, ligand_traj.top.residue(0).name)
            oechem.OEAddMols(receptor, ligand_posed)
            print("reset")

        else:
            ligand_posed = lig_initial
            print("stable")

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
        complex_prot = __eraseLigand(complex_prot, ligand_traj.top.residue(0).name)
        i = i + 1
    # Produce output trajectories for all combos of ligands
    output_dummy = tempfile.NamedTemporaryFile(suffix=".pdb")
    outputter = oechem.oemolostream(output_dummy.name)
    oechem.OEWriteMolecule(outputter, complex_prot)
    for ligand_combo in itertools.product(*ligand_stream):
        new_frame_prot = mdtraj.load(output_dummy.name)
        for built_ligand in ligand_combo:
            # Rename ligand
            ligand_dummy = tempfile.NamedTemporaryFile(suffix=".pdb")
            outputter_lig = oechem.oemolostream(ligand_dummy.name)
            oechem.OEWriteMolecule(outputter_lig, built_ligand)
            ligand_final_traj = mdtraj.load(ligand_dummy.name)
            for residue in ligand_final_traj.top.residues:
                if residue.name == 'UNL':
                    residue.name = built_ligand.GetTitle()
                    print(residue.name)
            new_frame_prot = new_frame_prot.stack(ligand_final_traj)
        output.append(new_frame_prot)
    print(output)
    print("End")
    return output


def __docking_internal(receptor: oechem.OEGraphMol, bound_ligand: oechem.OEGraphMol, docking_ligand: oechem.OEGraphMol):
    """
    """
    # Create posit config
    oedocking.OEReceptorSetBoundLigand(receptor, bound_ligand)
    poser = oedocking.OEPosit()
    poser.Initialize(receptor)

    # Create multiple conformations
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfs(100)
    omega_driver = oeomega.OEOmega(omegaOpts)
    conformer_docking = oechem.OEMol()  # type: OEMol
    print(oechem.OESmilesToMol(conformer_docking, NameToSMILES.addH(oechem.OEMolToSmiles(docking_ligand))))
    if not omega_driver(conformer_docking):
        single_sdf = tempfile.NamedTemporaryFile(suffix=".sdf")
        double_sdf = tempfile.NamedTemporaryFile(suffix=".sdf")
        smiles = NameToSMILES.addH(oechem.OEMolToSmiles(docking_ligand))
        print(smiles)
        subprocess.run("obabel -:'%s' -O %s --gen3d; obabel %s -O %s --confab --conf = 100" % (
        smiles, single_sdf.name, single_sdf.name, double_sdf.name), shell=True)
        ins = oechem.oemolistream()
        ins.SetConfTest(oechem.OEOmegaConfTest())
        ins.open(double_sdf.name)
        print(oechem.OEReadMolecule(ins, conformer_docking))

    # Failed conformer generation
    if conformer_docking.NumConfs() <= 1:
        raise TypeError('Conformer generation failed. Only %s conformers present' % conformer_docking.NumConfs())
    # Dock and get top conformer
    posed_ligand = oechem.OEGraphMol()  # type: OEGraphMol
    poser.DockMultiConformerMolecule(posed_ligand, conformer_docking)
    return posed_ligand


def __eraseLigand(complex, ligandCode):
    for atom in complex.GetAtoms():
        if oechem.OEAtomGetResidue(atom).GetName() == ligandCode:
            complex.DeleteAtom(atom)
    return complex
