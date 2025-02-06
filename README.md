# GROMACS Simulation Workflow for Protein-Ligand Complex

This repository provides a detailed workflow for performing molecular dynamics (MD) simulations using GROMACS, including the preparation of protein-ligand complexes and the steps required for running simulations. The workflow is inspired by GROMACS tutorials by Justin A. Lemkul, Ph.D., from the Virginia Tech Department of Biochemistry.

## Workflow Overview

The steps in this workflow guide the user through the process of preparing and running molecular dynamics simulations of protein-ligand complexes, starting from structure preparation to energy analysis.

---

## Commands:

```bash
# Protein structure preparation

# Convert PDB to GRO file
$ gmx pdb2gmx -f protein_clean.pdb -o protein_processed.gro -ter

# Fix ligand bonds
$ perl sort_mol2_bonds.pl lig.mol2 lig_fix.mol2

# Prepare ligand using CHARMM force field
$ python cgenff_charmm2gmx.py JZ4 lig_fix.mol2 jz4.str charmm36-jul2022.ff

# Process and edit ligand/protein structures
$ gmx editconf -f lig_ini.pdb -o lig.gro
$ gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0

# Add solvent molecules
$ gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

# Ionization
$ gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
$ gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Energy minimization
$ gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
$ gmx mdrun -v -deffnm em

# Create index file for ligand
$ gmx make_ndx -f lig.gro -o index_ligand.ndx
> 0 & ! a H*
$ gmx genrestr -f lig.gro -n index_jz4.ndx -o posre_lig.itp -fc 1000 1000 1000

# Create general index file
$ gmx make_ndx -f em.gro -o index.ndx
> 1 | 13

# Equilibration: NVT
$ gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
$ gmx mdrun -deffnm nvt

# Equilibration: NPT
$ gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
$ gmx mdrun -deffnm npt

# Production run
$ gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
$ gmx mdrun -deffnm md_0_10

# Trajectory processing
$ gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
$ gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o start.pdb -dump 0
$ gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o md_0_10_fit.xtc -fit rot+trans

# Distance calculation
$ gmx distance -s md_0_10.tpr -f md_0_10_center.xtc -select 'resname "ligand" and name S plus resid 105 and name NE1' -oall

# Interaction groups using PyMOL and index files
$ gmx make_ndx -f em.gro -o index.ndx
> 13 & a S | a NE1
> 1 & r 105 & a HE1
> 17 | 18

# RMSD analysis
$ gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd_ligand.xvg

# Interaction energy calculation
$ gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr
$ gmx mdrun -deffnm ie -rerun md_0_10.xtc -nb cpu
$ gmx energy -f ie.edr -o interaction_energy.xvg
$ gmx hbond -s md_0_10.tpr -f md_0_10.xtc -n index.ndx -num hbonds.xvg

# Cluster analysis
$ gmx_mpi cluster -f md_0_10.xtc -s md_0_10.tpr -o cluster.xpm -g cluster.log -n index.ndx -cutoff 0.15 -method gromos -cl rep_structures.pdb 
$ echo 1 | gmx_mpi trjconv -s md_0_10.tpr -f md_0_10.xtc -o largest_cluster.pdb -n index.ndx

# Covariance Matrix Calculation
$ gmx covar -s md_0_10.tpr -f md_0_10_center.xtc -n index.ndx -o eigenvalues.xvg -v eigenvectors.trr -av average.pdb


## Steps

### 1. Protein Structure Preparation:
- Convert the protein structure from a PDB file to a GRO file using `pdb2gmx`.
- Process and fix ligand files using custom Perl and Python scripts.
- Prepare the complex for simulation by editing the protein-ligand structure using `editconf`.

### 2. Box and Solvent Addition:
- Define the simulation box with a dodecahedron shape and 1.0 nm distance from the solute using `editconf`.
- Add water molecules to the simulation box using `solvate`.

### 3. Ionization:
- Generate a topology file and prepare the system for ionization.
- Perform ionization using `genion`, ensuring the system is neutralized.

### 4. Energy Minimization (EM):
- Perform energy minimization to relax the system using `grompp` and `mdrun`.

### 5. Equilibration:
- Perform temperature (NVT) and pressure (NPT) equilibration using the appropriate `.mdp` files.
- Use `grompp` to prepare the systems and `mdrun` for simulation.

### 6. Production Run:
- Run the production molecular dynamics simulation (md), then process the trajectory using `trjconv` to analyze the system’s behavior.

### 7. Analysis:
- Perform various analyses including RMSD, distance calculations, and interaction energy evaluations using tools like `rms`, `distance`, and `energy`.

### 8. Hydrogen Bond Analysis
-Calculate the hydrogen bonds formed during the simulation and output the results in the hbonds.xvg file, which contains the number of hydrogen bonds over time.
---

## 9. Pressure Analysis

### Compute Pressure Autocorrelation
```bash
$ gmx_mpi analyze -f pressuren.xvg -ac
```

### Plot Pressure Data Using Gnuplot
```bash
gnuplot -e "set terminal png size 800,600; set output 'pressure_with_reference.png'; set title 'Pressure vs Time'; set xlabel 'Time (ps)'; set ylabel 'Pressure (bar)'; set grid; plot 'pressuren.xvg' using 1:2 with lines title 'Pressure' linecolor rgb 'purple', 'pressuren.xvg' using 1:2 smooth sbezier with lines title '10-ps Running Avg' linecolor rgb 'red';"
```

---

## Key Commands Used:
- **`gmx pdb2gmx`**: Converts PDB to GRO for GROMACS compatibility.
- **`perl sort_mol2_bonds.pl`**: Fixes the ligand structure.
- **`python cgenff_charmm2gmx.py`**: Prepares the ligand using the CHARMM force field.
- **`gmx solvate`**: Adds solvent molecules to the simulation box.
- **`gmx grompp`**: Prepares the system for energy minimization, equilibration, and production runs.
- **`gmx mdrun`**: Runs the simulation.
- **`gmx trjconv`**: Processes the trajectory for analysis.
- **`gmx rms`, `gmx distance`, `gmx energy`**: Perform analysis on the simulation data.

---

## File Descriptions:
- **`protein_clean.pdb`**: Cleaned protein structure.
- **`ligand.mol2`**: Ligand in MOL2 format.
- **`ligand_fix.mol2`**: Processed ligand structure.
- **`topol.top`**: Topology file for the system.
- **`em.mdp`, `nvt.mdp`, `npt.mdp`**: Parameter files for energy minimization, NVT, and NPT simulations.
- **`md.mdp`**: Parameter file for the production MD run.

---

## Interaction Energy Analysis:

After obtaining the `interaction_energy.xvg` file from GROMACS, you can use the `xvg.py` Python script for further analysis. The script processes the interaction energy data and computes energy components like Coulombic and Van der Waals interactions over time.

### Example Command:
```bash
$ python xvg.py 
```

The script outputs the interaction energy data in a format suitable for further analysis or plotting.

---

## Prerequisites:
- **GROMACS 2018** or later
- **Python 3.x**
- **Perl**
- **CHARMM36 force field files**

Ensure GROMACS is installed and set up correctly. Install Python and required libraries (e.g., MDAnalysis). Ensure Perl is installed for running the custom scripts.

---

## Conclusion:

This repository provides a systematic approach for simulating protein-ligand interactions using GROMACS, offering insights into molecular dynamics simulations and analysis. The workflow is modular, allowing you to adapt the steps for your specific system or research.

