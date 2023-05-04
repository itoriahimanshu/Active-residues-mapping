import parmed as pmd
from simtk import openmm, unit
from simtk.openmm import app
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the PDB file
pdb = app.PDBFile('myprotein.pdb')

# Create a ParmEd structure from the PDB file
structure = pmd.load_file('myprotein.pdb')

# Identify the active residues (in this example, residues with IDs 38 and 46)
active_residues = [38, 46]

# Create an OpenMM system object using an MM force field
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=openmm.app.NoCutoff)

# Compute the partial charges of each atom in the protein
integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
platform = openmm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
openmm.LocalEnergyMinimizer.minimize(simulation.context)

# Compute the total charge of each active residue
charges = []
for residue in structure.residues:
    if residue.number in active_residues:
        charge = 0.0
        for atom in residue.atoms:
            force = system.getForce(atom.forceGroup)
            if isinstance(force, openmm.NonbondedForce):
                charge += atom.charge.value_in_unit(unit.elementary_charge)
        charges.append(charge)

# Get the coordinates of the active residues
coordinates = []
for atom in structure.atoms:
    if atom.residue.number in active_residues:
        x = atom.xx
        y = atom.xy
        z = atom.xz
        coordinates.append([x, y, z])

# Create a 3D scatter plot with the charges as colors
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.array(coordinates)[:, 0], np.array(coordinates)[:, 1], np.array(coordinates)[:, 2], c=charges)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
