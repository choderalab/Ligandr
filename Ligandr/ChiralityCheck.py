"""Partial clone of the VMD plug-ins `cispeptide` and `chirality` """
""" 
Directly taken from Josh Fass- find the original at 
https://github.com/choderalab/msm-pipeline/blob/utils/msmpipeline/utils.py
"""

import mdtraj as md
import numpy as np


def find_peptide_bonds(top):
    """Find all peptide bonds in a topology
    Parameters
    ----------
    top : mdtraj topology
    Returns
    -------
    peptide_bonds : list of 4-tuples
        each 4-tuple contains indices of the (CA, N, C, O) atoms involved in the peptide bond
    residues : list of integers
        peptide_bonds[i] involves resids (residues[i], residues[i]+1)
    References
    ----------
    [1] Stereochemical errors and their implications for molecular dynamics simulations.
        Schriener et al., 2011. BMC Bioinformatics. DOI: 10.1186/1471-2105-12-190
        http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-190
    """

    bond_indices = [sorted((bond[0].index, bond[1].index)) for bond in top.bonds]

    peptide_bonds = []
    residues = []

    for i in range(top.n_residues - 1):

        j = i + 1

        # each of these selections returns an array containing zero or one integer
        O_ind = (top.select('resid=={0} and name == O'.format(i)))
        C_ind = (top.select('resid=={0} and name == C'.format(i)))

        N_ind = (top.select('resid=={0} and name == N'.format(j)))
        CA_ind = (top.select('resid=={0} and name == CA'.format(j)))

        # if C(i) is bonded to N(j), and all of the index arrays are of length-1, add the (CA, N, C, O) tuple
        # to our peptide_bonds list
        if sorted((C_ind, N_ind)) in bond_indices and sum([len(a) != 1 for a in (O_ind, C_ind, N_ind, CA_ind)]) == 0:
            peptide_bonds.append((CA_ind[0], N_ind[0], C_ind[0], O_ind[0]))
            residues.append(i)

    return peptide_bonds, residues


def check_cispeptide_bond(traj, dihedral_cutoff=85):
    """Given a trajectory, check every peptide bond in every frame.
    Return a boolean array of shape=(len(traj), n_peptide_bonds)
    Print a warning if any cis-peptide bonds are found
    Parameters
    ----------
    traj : mdtraj.Trajectory
    dihedral_cutoff : default 85
        Dihedral angle cutoff, in degrees
    Returns
    -------
    problems : boolean array
        problems[i,j] means that, in frame i, peptide bond j is cis
    References
    ----------
    [1] Stereochemical errors and their implications for molecular dynamics simulations.
        Schriener et al., 2011. BMC Bioinformatics. DOI: 10.1186/1471-2105-12-190
        http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-190
    """

    indices, residues = find_peptide_bonds(traj.top)
    angles = md.compute_dihedrals(traj, indices)
    in_degrees = angles * 180 / np.pi
    problems = in_degrees > dihedral_cutoff

    if np.sum(problems) > 0:
        print('Problems in {0} frames!'.format(sum(problems.sum(1) != 0)))
        print('Problems in {0} bonds!'.format(sum(problems.sum(0) != 0)))

    return problems


def define_chiral_centers():
    """Constructs a dictionary of the chiral centers of each amino acid
    Returns
    -------
    chiral_dict : dictionary
        maps three-letter residue names to lists of 4-tuples of atom *names*
    Notes / to-do's
    ---------------
    - to-do: find an independent reference here (modified this from
      the definition of chiral centers from lines 53-63 in chirality.tcl in VMD)
    - to-do: check correctness
    References
    ----------
    [1] Stereochemical errors and their implications for molecular dynamics simulations.
        Schriener et al., 2011. BMC Bioinformatics. DOI: 10.1186/1471-2105-12-190
        http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-190
    """
    chiral_dict = dict()
    AAs = 'ALA ARG ASN ASP CYS GLN GLU HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()
    for residue_name in AAs:
        chiral_dict[residue_name] = [tuple('HA N C CB'.split()),
                                     tuple('CA N C CB'.split())]

        if residue_name in 'THR ILE'.split():
            chiral_dict[residue_name].append(tuple('HB CA OG1 CG2'.split()))
            chiral_dict[residue_name].append(tuple('CB CA OG1 CG2'.split()))

    return chiral_dict


def get_chiral_atoms(top):
    """Given a topology, return a list of 4-tuples of chiral atom indices
    We'll then compute dihedral angles using each 4-tuple, and check if any are below 0
    Parameters
    ----------
    top : mdtraj.topology
    Returns
    -------
    chiral_atoms : list of 4-tuples
    """

    # dictionary mapping from 3-letter residue name to list of 4-tuples of atom *names*
    chiral_dict = define_chiral_centers()

    # a list containing 4-tuples of atom *indices* within top, defining torsions we'll check
    chiral_atoms = []

    for residue in top.residues:
        # if this is a residue that contains chiral positions
        if residue.name in chiral_dict.keys():

            # get all the atoms in the current residue
            atoms = [a for a in residue.atoms]
            atom_names = [a.name for a in atoms]

            # for each quartet, retrieve the corresponding atom indices within top
            for name_quartet in chiral_dict[residue.name]:
                index_quartet = []
                for name in name_quartet:
                    current_atom = atoms[atom_names.index(name)]
                    index_quartet.append(current_atom.index)
                chiral_atoms.append(tuple(index_quartet))

    return chiral_atoms


def check_chirality(traj, threshold=0):
    """Check whether each chiral amino acid is in the L-configuration
    at each frame in the trajectory.
    Parameters
    ----------
    traj : mdtraj.trajectory
    threshold : angle (in degrees)
        angles below this threshold will be taken as errors
    Returns
    -------
    chiral_atoms : list of 4-tuples
    errors : numpy.ndarray, shape=(len(traj), len(chiral_atoms)), dtype=bool
        errors[i,j] means there is a chirality error in frame i, quartet j
    References
    ----------
    [1] Stereochemical errors and their implications for molecular dynamics simulations.
        Schriener et al., 2011. BMC Bioinformatics. DOI: 10.1186/1471-2105-12-190
        http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-190
    """
    chiral_atoms = get_chiral_atoms(traj.top)
    dihedrals = md.compute_dihedrals(traj, chiral_atoms)
    print(dihedrals)
    in_degrees = dihedrals * 180 / np.pi
    errors = in_degrees < threshold
    return chiral_atoms, errors
