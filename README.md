# GROMACS Simulation Workflow for Protein-Ligand Complex

This repository provides a detailed workflow for performing molecular dynamics (MD) simulations using GROMACS, including the preparation of protein-ligand complexes and the steps required for running simulations. The workflow is inspired by GROMACS tutorials by **Justin A. Lemkul, Ph.D.**, from the Virginia Tech Department of Biochemistry.

---

## Workflow Overview

The steps in this workflow guide the user through the process of preparing and running molecular dynamics simulations of protein-ligand complexes, starting from structure preparation to energy analysis.

---

## Commands

### Protein Structure Preparation
1. Convert the protein structure from a PDB file to a GRO file:
   ```bash
   gmx pdb2gmx -f protein_clean.pdb -o protein_processed.gro -ter

Process and fix ligand files using custom Perl and Python scripts:
perl sort_mol2_bonds.pl lig.mol2 lig_fix.mol2
python cgenff_charmm2gmx.py JZ4 lig_fix.mol2 jz4.str charmm36-jul2022.ff

Edit the complex structure:
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

Box and Solvent Addition
Add solvent molecules:

gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
Ionization
Neutralize the system:

gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
Energy Minimization (EM)
Minimize the energy:

gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
Equilibration
Perform NVT and NPT equilibration:

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
Production Run
Run the production simulation:

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
gmx mdrun -deffnm md_0_10
Process the trajectory:

gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
Analysis
RMSD Analysis:

gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd_ligand.xvg
Interaction Energy Analysis:

gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr
gmx mdrun -deffnm ie -rerun md_0_10.xtc -nb cpu
gmx energy -f ie.edr -o interaction_energy.xvg

Visualization with Gnuplot
After analyzing interaction energies, you can visualize the results using Gnuplot.

Example: Plotting Interaction Energies
Run the following Gnuplot command to generate a plot:

gnuplot -e "
set terminal png size 800,600;
set output 'interactionenergy_with_reference.png';
set title 'Energy vs Time';
set xlabel 'Time (ps)';
set ylabel 'Energy (kJ/mol)';
set grid;
plot 'interaction_energy.xvg' using 1:2 with lines title 'Pressure' linecolor rgb 'purple', \
     'interaction_energy.xvg' using 1:2 smooth sbezier with lines title '10-ps Running Avg' linecolor rgb 'red';"
This will create a file named interactionenergy_with_reference.png with the following details:
X-Axis: Time (ps)
Y-Axis: Energy (kJ/mol)
Purple Line: Raw data.
Red Line: 10-ps running average (smoothed).

File Descriptions
protein_clean.pdb: Cleaned protein structure.
ligand.mol2: Ligand in MOL2 format.
ligand_fix.mol2: Processed ligand structure.
topol.top: Topology file for the system.
em.mdp, nvt.mdp, npt.mdp: Parameter files for energy minimization, NVT, and NPT simulations.
md.mdp: Parameter file for the production MD run.

Prerequisites
GROMACS 2018 or later
Python 3.x
Perl
CHARMM36 force field files
Gnuplot for visualization
Make sure all dependencies are installed and your environment is properly configured before starting.

Conclusion
This repository provides a systematic approach for simulating protein-ligand interactions using GROMACS, offering insights into molecular dynamics simulations and analysis. The workflow is modular, allowing you to adapt the steps for your specific system or research.
