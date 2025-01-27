GROMACS Simulation Workflow for Protein-Ligand Complex
This repository provides a detailed workflow for performing molecular dynamics (MD) simulations using GROMACS, including the preparation of protein-ligand complexes and the steps required for running simulations. The workflow is inspired by GROMACS tutorials by Justin A. Lemkul, Ph.D., from Virginia Tech Department of Biochemistry.

Workflow Overview
The steps in this workflow guide the user through the process of preparing and running molecular dynamics simulations of protein-ligand complexes, starting from structure preparation to energy analysis.

Steps:
Protein Structure Preparation:

Convert the protein structure from a PDB file to a GRO file using pdb2gmx.
Process and fix ligand files using custom Perl and Python scripts.
Prepare the complex for simulation by editing the protein-ligand structure using editconf.
Box and Solvent Addition:

Define the simulation box with a dodecahedron shape and 1.0 nm distance from the solute using editconf.
Add water molecules to the simulation box using solvate.
Ionization:

Generate a topology file and prepare the system for ionization.
Perform ionization using genion, ensuring the system is neutralized.
Energy Minimization (EM):

Perform energy minimization to relax the system using grompp and mdrun.
Equilibration:

Perform temperature (NVT) and pressure (NPT) equilibration using the appropriate .mdp files.
Use grompp to prepare the systems and mdrun for simulation.
Production Run:

Run the production molecular dynamics simulation (md), then process the trajectory using trjconv to analyze the system’s behavior.
Analysis:

Perform various analyses including RMSD, distance calculations, and interaction energy evaluations using tools like rms, distance, and energy.

Pressure Analysis:
Use gmx_mpi analyze to compute the autocorrelation function for pressure.
gmx_mpi analyze -f pressuren.xvg -ac

Plot the pressure data using gnuplot with a 10-ps running average.

gnuplot -e "
set terminal png size 800,600;
set output 'pressure_with_reference.png';
set title 'Pressure vs Time';
set xlabel 'Time (ps)';
set ylabel 'Pressure (bar)';
set grid;
plot 'pressuren.xvg' using 1:2 with lines title 'Pressure' linecolor rgb 'purple', \
     'pressuren.xvg' using 1:2 smooth sbezier with lines title '10-ps Running Avg' linecolor rgb 'red';"

Key Commands Used:
gmx pdb2gmx: Converts PDB to GRO for GROMACS compatibility.
perl sort_mol2_bonds.pl: Fixes the ligand structure.
python cgenff_charmm2gmx.py: Prepares the ligand using the CHARMM force field.
gmx solvate: Adds solvent molecules to the simulation box.
gmx grompp: Prepares the system for energy minimization, equilibration, and production runs.
gmx mdrun: Runs the simulation.
gmx trjconv: Processes the trajectory for analysis.
gmx rms, gmx distance, gmx energy: Perform analysis on the simulation data.

File Descriptions:
protein_clean.pdb: Cleaned protein structure.
ligand.mol2: Ligand in MOL2 format.
ligand_fix.mol2: Processed ligand structure.
topol.top: Topology file for the system.
em.mdp, nvt.mdp, npt.mdp: Parameter files for energy minimization, NVT, and NPT simulations.
md.mdp: Parameter file for the production MD run.

Interaction Energy Analysis:

After obtaining the interaction_energy.xvg file from GROMACS, you can use the xvg.py Python script for further analysis.
The xvg.py script processes the interaction energy data from the .xvg file and computes various energy components, such as Coulombic and Van der Waals interactions, over time.
Example command to use the xvg.py script:

python xvg.py interaction_energy.xvg
The script will output the interaction energy data in a format suitable for further analysis or plotting.

Prerequisites:
GROMACS 2018 or later
Python 3.x
Perl
CHARMM36 force field files


Ensure GROMACS is installed and set up correctly.
Install Python and required libraries (e.g., MDAnalysis or any other Python dependencies you may need).
Ensure Perl is installed for running the custom scripts.
Follow the steps in the tutorial to prepare and run your simulation.

Conclusion:
This repository provides a systematic approach for simulating protein-ligand interactions using GROMACS, offering insights into molecular dynamics simulations and analysis. The workflow is modular, allowing you to adapt the steps for your specific system or research.
