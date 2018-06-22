"""
Overlays two protein chains and transfers ligands and other bound molecule
"""

import mdtraj
from Bio import pairwise2
import numpy
from typing import Dict, Tuple, List, Union
import itertools
from copy import deepcopy


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

    # Find steric clashes
    ligand_residues = [residue.index for residue in ligand_added_target.atom_slice(
        ligand_added_target.top.select(ligand_code)).top.residues]  # Find the residue of the ligand
    other_residues = [residue.index for residue in posed_target.top.residues]
    pairs = list(itertools.product(other_residues, ligand_residues))
    distances = mdtraj.compute_contacts(ligand_added_target, pairs, scheme='closest', ignore_nonprotein=False)[0]
    steric_clash = []
    for index, frame in enumerate(distances):
        min_distance = numpy.amin(frame)
        if min_distance * 10 <= threshold:
            steric_clash.append((index, min_distance * 10))  # Convert to angstroms
    print(steric_clash)
    return ligand_added_target, steric_clash
