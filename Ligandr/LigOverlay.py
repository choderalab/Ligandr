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
from openeye import oechem, oedocking, oeomega
from openeye.oechem import OEMol, OEGraphMol
import itertools
import subprocess
import parmed
from io import StringIO
import os.path


def dock_molecule(trajectory_frame, original_ligands, new_ligands, title_list, keep_ions=None,
                  override_replacement=False, gen_xml=False, location="./") -> List[
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

    gen_xml: bool, optional, default:False
        Setting to true runs the generate xml subroutine and saves openmm forcefield compatible XMLs in location. WARNING:
        Using this option requires antechmaber and is time-consuming

    location: str, optional, default:None
        The location to save xml files. Default is current directory
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
    complex_prot = _trajectory_to_oemol_(overall_dry_frame)

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
        # Find approximate location of binding pocket
        ligand_center = np.average(ligand_traj.xyz[0], axis=0) * 10  # Convert from nm to angstrom

        # Explicit value necessary for generated methods
        x = float(ligand_center[0])
        y = float(ligand_center[1])
        z = float(ligand_center[2])

        # Create receptor molecule
        receptor = oechem.OEGraphMol()
        oedocking.OEMakeReceptor(receptor, complex_prot, x, y, z)

        # Write out and open ligand in OpenEye
        lig_initial = _trajectory_to_oemol_(ligand_traj)
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
                ligand_posed = _docking_internal(receptor, lig_initial, validation_ligand)
                complex_prot = _trajectory_to_oemol_(
                    overall_dry_frame.atom_slice(overall_dry_frame.top.select("not ( %s )" % ligand)))
                oechem.OEAddMols(receptor, ligand_posed)
                print(
                    "Warning: Original ligand had either incorrect stereochemistry or missing heavy atoms. Ligand from pdb database placed instead. ")

            else:
                ligand_posed = lig_initial
        else:
            print("ligand not found in PDB: No validation and continuing")
        # Dock new ligands
        ligand_position = []  # type: List[OEGraphMol]
        j = 0
        for new_ligand in new_ligands[i]:
            new_structure = oechem.OEGraphMol()
            oechem.OESmilesToMol(new_structure, new_ligand)
            new_structure.SetTitle("ligand replacement")
            print(title_list[i][j])
            docked_ligand = _docking_internal(receptor, ligand_posed, new_ligand)
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
            ligand_named = tempfile.NamedTemporaryFile(suffix=".pdb")
            ligand_final_traj.save(ligand_named.name)
            if gen_xml:
                if not os.path.isfile(location + built_ligand.GetTitle() + ".xml"):
                    generate_xml(ligand_named.name, location, oechem.OENetCharge(built_ligand), built_ligand.GetTitle())
                else:
                    print(built_ligand.GetTitle() + ".xml already exists")
            new_frame_prot = new_frame_prot.stack(ligand_final_traj)
        output.append(new_frame_prot)
    print("End")
    return output


def generate_xml(ligand_pdb, location, charge, title="UNK"):
    file_stub = location + title
    subprocess.run("antechamber -i %s -fi pdb -o %s.mol2 -fo mol2 -c bcc -nc %s" % (ligand_pdb, file_stub, charge),
                   shell=True)
    subprocess.run("parmchk2 -i %s.mol2 -f mol2 -o %s.frcmod" % (file_stub, file_stub), shell=True)
    params_mol = parmed.amber.AmberParameterSet.from_leaprc(
        StringIO(u'SAM = loadMol2 %s.mol2\nloadAmberParams %s.frcmod' % (file_stub, file_stub)))
    params_sam = parmed.openmm.OpenMMParameterSet.from_parameterset(params_mol)
    params_sam.write('%s.xml' % file_stub)
    print("Forcefield generation of %s complete" % title)


def _docking_internal(receptor: oechem.OEGraphMol, bound_ligand: oechem.OEGraphMol, ligand_string):
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
    An oechem molecule with coordinates of its docking

    """
    # Create posit config
    oedocking.OEReceptorSetBoundLigand(receptor, bound_ligand)
    poser = oedocking.OEHybrid(oedocking.OEDockMethod_Hybrid2, oedocking.OESearchResolution_High)
    poser.Initialize(receptor)

    # Create multiple conformations
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfs(1000)
    omega_driver = oeomega.OEOmega(omegaOpts)
    conformer_docking = oechem.OEMol()  # type: OEMol
    oechem.OESmilesToMol(conformer_docking, ligand_string)
    oechem.OEAddExplicitHydrogens(conformer_docking)
    print(conformer_docking.NumAtoms())

    if not omega_driver(conformer_docking):
        single_sdf = tempfile.NamedTemporaryFile(suffix=".sdf")
        double_sdf = tempfile.NamedTemporaryFile(suffix=".sdf")
        smiles = ligand_string
        subprocess.run("obabel -:'%s' -O %s --gen3d; obabel %s -O %s --confab --conf = 100" % (
            smiles, single_sdf.name, single_sdf.name, double_sdf.name), shell=True)
        ins = oechem.oemolistream()
        ins.SetConfTest(oechem.OEOmegaConfTest())
        ins.open(double_sdf.name)
        oechem.OEReadMolecule(ins, conformer_docking)

    # Dock and get top conformer
    posed_ligand = oechem.OEGraphMol()  # type: OEGraphMol
    poser.DockMultiConformerMolecule(posed_ligand, conformer_docking)
    return posed_ligand


def _trajectory_to_oemol_(trajectory):
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
