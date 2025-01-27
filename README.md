GROMACS Simulation Workflow for Protein-Ligand Complex
This repository provides a detailed workflow for performing molecular dynamics (MD) simulations using GROMACS, including the preparation of protein-ligand complexes and the steps required for running simulations. The workflow is inspired by GROMACS tutorials by Justin A. Lemkul, Ph.D., from Virginia Tech Department of Biochemistry.

Workflow Overview
The steps in this workflow guide the user through the process of preparing and running molecular dynamics simulations of protein-ligand complexes, starting from structure preparation to energy analysis.

Commands:
gmx pdb2gmx -f ptotein_clean.pdb -o protein_processed.gro -ter
perl sort_mol2_bonds.pl lig.mol2 lig_fix.mol2
python cgenff_charmm2gmx.py JZ4 lig_fix.mol2 jz4.str charmm36-jul2022.ff
gmx editconf -f lig_ini.pdb -o lig.gro
gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
gmx make_ndx -f lig.gro -o index_ligand.ndx
...
 > 0 & ! a H*
 > q

gmx genrestr -f lig.gro -n index_jz4.ndx -o posre_lig.itp -fc 1000 1000 1000
gmx make_ndx -f em.gro -o index.ndx
> 1 | 13
> q

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr

gmx mdrun -deffnm nvt

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr

gmx mdrun -deffnm npt

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr

gmx mdrun -deffnm md_0_10

gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact

gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o start.pdb -dump 0

gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o md_0_10_fit.xtc -fit rot+trans

gmx distance -s md_0_10.tpr -f md_0_10_center.xtc -select 'resname "ligand" and name S plus resid 105 and name NE1' -oall(These atoms were determined based on the residues that form interactions in the binding region between the ligand and the protein, as well as the ligand atoms by using Pymol).

gmx make_ndx -f em.gro -o index.ndx
 > 13 & a S | a NE1
 (creates group 17)
 > 1 & r 105 & a HE1
 (creates group 18)
 >17 | 18
 > q

gmx make_ndx -f em.gro -n index.ndx
...
 > 13 & ! a H*
name 20 ligand_Heavy

gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd_ligand.xvg

gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr

gmx mdrun -deffnm ie -rerun md_0_10.xtc -nb cpu

gmx energy -f ie.edr -o interaction_energy.xvg 
GROMACS Tutorials
Justin A. Lemkul, Ph.D.
Virginia Tech Department of Biochemistry

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
