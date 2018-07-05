"""
Overlays two protein chains and transfers ligands and other bound molecule
"""

import mdtraj
from Bio import pairwise2
import numpy
from typing import Dict, Tuple, List, Union
from copy import deepcopy
import tempfile
from simtk import openmm, unit
from simtk.openmm import app

import openmmtools
from openmmtools import integrators, states, mcmc
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion, AlchemicalState


def align_structs(target_trajectory: mdtraj.Trajectory, reference_structure: mdtraj.Trajectory,
                  subset_target: str = "protein",
                  subset_reference: str = "protein",
                  atom_selection: str = "name CA") -> mdtraj.Trajectory:
    """
    Poses a target protein to a reference structure. The reference is based on the COORDINATES OF THE FIRST CHAIN for each
    structure.

    Parameters
    ----------
    reference_structure: mdtraj.Trajectory, required
        The reference structure to pose to. If this has multiple frames, the reference will be taken from frame 0.
    subset_target: str, optional, default = "protein"
        The string representing which residues to select for the target. This should be mdtraj compatible, and if this is
        filled out, then subset_reference will in general be filled up as well.
    subset_reference: str, optional, default = "protein"
        The string representing which residues to select for the reference. This should be mdtraj compatible, and if this is
        filled out, then subset_target will in general be filled up as well.
    atom_selection: str, optional default = "name CA"
        The type of atom to use to compare. It should be mdtraj compatible, and it is highly recommended to use the default.
    target_trajectory: mdtraj.Trajectory, required
        The target trajectory that will be transformed.

    Returns
    -------
    mdtraj.Trajectory
        A Trajectory object with coordinates transformed to the reference structure frame
    """
    s1 = target_trajectory.top.select("%s and %s" % (subset_target, atom_selection))
    s2 = reference_structure.top.select("%s and %s" % (subset_reference, atom_selection))
    fast_s1 = set(s1)
    fast_s2 = set(s2)
    # If the lengths of the original selection is the same, an alignment will probably work (Don't waste computation)
    if len(s1) != len(s2):
        resvalue_1, resvalue_2 = _find_alignment(target_trajectory.top.to_fasta()[0],
                                                 reference_structure.top.to_fasta()[
                                                     0])  # type: Tuple[List[int], List[int]]
        resvalue_1 = set(resvalue_1)
        resvalue_2 = set(resvalue_2)
        index_target = numpy.asarray([atom.index for atom in target_trajectory.top.atoms if
                                      (atom.residue.index in resvalue_1 and atom.index in fast_s1)])
        index_ref = numpy.asarray([atom.index for atom in reference_structure.top.atoms if
                                   (atom.residue.index in resvalue_2 and atom.index in fast_s2)])
        print(len(index_ref))
        print(len(index_target))
        assert len(index_ref) == len(index_target), "Different number of atoms selected from alignment"
        s1 = index_target
        s2 = index_ref
    print(len(s1))
    print(len(s2))
    overlay_targ = target_trajectory.superpose(reference_structure, atom_indices=s1, ref_atom_indices=s2)
    return overlay_targ


def _find_alignment(fasta_seq1: str, fasta_seq2: str) -> Union[List[int], List[int]]:
    """
    Given two fasta sequences, align and return the resid range for each matching area. End gaps are ignored for scoring.

    Parameters
    ----------
    fasta_seq1: str, required
        A string representing the amino acid code sequence of a protein

    fasta_seq2: str, required
        A string representing the amino acid code sequence of a protein. This will be aligned with the first string

    Returns
    -------
    List[int], List[int]
        Two lists of ints, representing the indices of the included residues for each sequence respectively
    """
    # Bookkeeping
    seq1_indices = []
    seq2_indices = []
    seq1_offset = 0
    seq2_offset = 0
    break1 = 0
    break2 = 0

    alignment = \
        pairwise2.align.globalxs(fasta_seq1, fasta_seq2, -0.5, -0.5, penalize_end_gaps=False, one_alignment_only=True)[
            0]
    assert len(alignment[0]) == len(alignment[1]), "Something went wrong with seq alignment."
    for i in range(0, len(alignment[0])):
        # Determine how long the opening sequence of gaps is
        if alignment[0][i] == '-':
            if not seq1_indices:
                seq1_offset += 1
            else:
                break1 += 1
        if alignment[1][i] == '-':
            if not seq2_indices:
                seq2_offset += 1
            else:
                break2 += 0

        # Find indices with a valid character for each object
        if alignment[0][i] != '-' and alignment[1][i] != '-':
            seq1_indices.append(i)
            seq2_indices.append(i)

    # Use the offset value
    if seq1_offset != 0:
        seq1_indices = [index - seq1_offset for index in seq1_indices]
    if seq2_offset != 0:
        seq2_indices = [index - seq2_offset for index in seq2_indices]
    return seq1_indices, seq2_indices


def place_ligand(posed_target: mdtraj.Trajectory, reference: mdtraj.Trajectory,
                 ligand_code="not protein and not water", threshold=2.0) -> Union[
    mdtraj.Trajectory, List[Tuple[int, int]]]:
    """
    Places a ligand and checks for steric clashes

    Parameters
    ----------
    threshold: double, optional, default = 2.0
        The threshold for steric clashes. The default is 2.0 angstroms.
    posed_target: mdtraj.Trajectory, required
        The target to place the ligand in. This must be longer than or equal to the number of frames in reference.
    reference: mdtraj.Trajectory, required
        The structure to extract the ligand from
    ligand_code: str, optional, default = "not protein and not water"
        An mdtraj selection string to select which molecules to place onto the target. Defaults to all non-solvent,
        non-protein molecules

    Returns
    -------

    """
    # Place ligand on trajectory
    assert len(reference.top.select(ligand_code)) > 0, "no ligand selected"
    if len(posed_target) != len(reference):
        dummy_reference = deepcopy(reference)
        for i in range(0, len(posed_target) - len(reference)):
            reference = reference.join(dummy_reference)
    ligand_added_target = posed_target.stack(reference.atom_slice(reference.top.select(ligand_code)))

    return ligand_added_target


def shrake_rupley_fractions(traj: mdtraj.Trajectory, ligand: str) -> numpy.array:
    """

    Parameters
    ----------
    traj: mdtraj.Trajectory, required
        The trajectory to measure the accessible area
    ligand: str, required
        A selection string specifying which ligand to measure accessible area for

    Returns
    -------
    numpy.array
        The total fraction of accessible area
    """
    just_ligand = traj.atom_slice(traj.top.select(ligand))  # type: mdtraj.Trajectory
    # Residue indexing gets messed up so save off
    traj_lig_dummy = tempfile.NamedTemporaryFile(suffix=".h5")
    just_ligand.save(traj_lig_dummy.name)
    just_ligand = mdtraj.load(traj_lig_dummy.name)
    traj_lig_dummy.close()

    # Calc shrake rupely
    ligand_surface = mdtraj.shrake_rupley(just_ligand, probe_radius=0, mode='residue')[:, 0]
    traj_surface_lig = mdtraj.shrake_rupley(traj, probe_radius=0, mode='residue')[:, -1]
    ligand_surface_frac = []
    for i in range(len(ligand_surface)):
        val = traj_surface_lig[i] / ligand_surface[i]
        ligand_surface_frac.append(val)

    ligand_surface_frac = numpy.asarray(ligand_surface_frac)
    return ligand_surface_frac


# Following code is heavily based off of JDC's iapetus


def alchemy_minimization(trajectory: mdtraj.Trajectory, ligand: str, forcefield_loader: List[str] = None,
                         step=1) -> Union(
    mdtraj.Trajectory, List[List[float]]):
    """

    Parameters
    ----------
    trajectory: mdtraj.Trajectory, required
        The trajectory to minimize
    ligand: str, required
        An mdtraj selection string that represents the area where to lower potential energy interactivity
    forcefield_loader: List[str], optional, default = None
        Various forcefields that can be loaded

    Returns
    -------
    mdtraj.Trajectory
        A trajectory where every frame has been minimized
    List[List[float]]
        A list composed of floats representing the energy at the following lambdas pre and post minimization: (0,0), (0.1,0)
        as well as at 1,1

    """

    # Bookkeeping
    output_traj = trajectory[0]
    total_energy = []

    # Create system
    trajectory = trajectory[::step]
    topology_main = trajectory.top.to_openmm()
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    for field in forcefield_loader:
        forcefield.loadFile(field)
    system = forcefield.createSystem(topology_main, nonbondedMethod=app.NoCutoff, rigidWater=True,
                                     ewaldErrorTolerance=0.0005)

    # Do something with LJ potentials
    for force in system.getForces():
        if force.__class__.__name__ == 'NonbondedForce':
            for index in range(system.getNumParticles()):
                [charge, sigma, epsilon] = force.getParticleParameters(index)
                if sigma / unit.nanometers == 0.0:
                    force.setParticleParameters(index, charge, 1.0 * unit.angstroms, epsilon)

    # Get alchemical system
    alchemy_system = _soften_ligand(trajectory.top.select(ligand), system)
    alchemical_state = AlchemicalState.from_system(alchemy_system)  # type: AlchemicalState
    thermo_state = states.ThermodynamicState(alchemy_system, temperature=300.0 * unit.kelvin)
    alchem_thermo_state = states.CompoundThermodynamicState(thermo_state, [alchemical_state])

    n_annealing_steps = 1000
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 1.0 * unit.femtoseconds)
    context, integrator = openmmtools.cache.global_context_cache.get_context(alchem_thermo_state,
                                                                             integrator)  # type: (openmm.Context, openmm.Integrator)
    # Run dynamics
    for frame in trajectory:
        energy = []
        sampler_state = states.SamplerState(positions=frame.xyz[0])
        sampler_state.apply_to_context(context)

        # Each of these steps calculate and store energy
        alchemical_state.lambda_sterics = 0.0
        alchemical_state.lambda_electrostatics = 0.0
        alchemical_state.apply_to_context(context)
        state_test = context.getState(getPositions=True, getEnergy=True)
        energy.append(state_test.getPotentialEnergy()._value)

        alchemical_state.lambda_sterics = 0.1
        alchemical_state.lambda_electrostatics = 0.0
        alchemical_state.apply_to_context(context)
        state_test = context.getState(getPositions=True, getEnergy=True)
        energy.append(state_test.getPotentialEnergy()._value)

        for step in range(n_annealing_steps):
            alchemical_state.lambda_sterics = float(step) / float(n_annealing_steps)
            alchemical_state.lambda_electrostatics = 0.0
            alchemical_state.apply_to_context(context)
            integrator.step(1)
        for step in range(n_annealing_steps):
            alchemical_state.lambda_sterics = 1.0
            alchemical_state.lambda_electrostatics = float(step) / float(n_annealing_steps)
            alchemical_state.apply_to_context(context)
            integrator.step(1)
        sampler_state.update_from_context(context)

        alchemical_state.lambda_sterics = 0.0
        alchemical_state.lambda_electrostatics = 0.0
        alchemical_state.apply_to_context(context)
        state_test = context.getState(getPositions=True, getEnergy=True)
        energy.append(state_test.getPotentialEnergy()._value)

        alchemical_state.lambda_sterics = 0.1
        alchemical_state.lambda_electrostatics = 0.0
        alchemical_state.apply_to_context(context)
        state_test = context.getState(getPositions=True, getEnergy=True)
        energy.append(state_test.getPotentialEnergy()._value)

        alchemical_state.lambda_sterics = 1.0
        alchemical_state.lambda_electrostatics = 1.0
        alchemical_state.apply_to_context(context)
        state_test = context.getState(getPositions=True, getEnergy=True)
        energy.append(state_test.getPotentialEnergy()._value)

        # Generate ouput frame
        positions = state_test.getPositions()
        fixed_positions = []
        fixed_positions.append([[coord._value for coord in xyz] for xyz in positions])
        traj_out = mdtraj.Trajectory(fixed_positions, topology=frame.top, unitcell_lengths=frame.unitcell_lengths[0],
                                     unitcell_angles=frame.unitcell_angles[0])
        output_traj = output_traj.join(traj_out)
        total_energy.append(energy)

    return output_traj[1:], total_energy


def _soften_ligand(ligand_index: List[int], reference_system: openmm.System) -> openmm.System:
    """

    Parameters
    ----------
    ligand_index: List[int], required
        The indecies of the atoms to soften
    reference_system: openmm.System, required
        The system that the ligand is part of

    Returns
    -------
    openmm.System
        A system ready for alchemical changes

    """
    factory = AbsoluteAlchemicalFactory(consistent_exceptions=False, disable_alchemical_dispersion_correction=True)
    alchemical_region = AlchemicalRegion(alchemical_atoms=ligand_index)
    alchemical_system = factory.create_alchemical_system(reference_system, alchemical_region)
    return alchemical_system
