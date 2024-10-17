## Overview
Mutational Energy Landscape Trap (MELT) is a computational tool designed to control protein dynamics by combining Normal Mode Analysis (NMA) and in silico mutagenesis. It offers a novel approach for manipulating protein structures, focusing on their dynamic behaviors rather than static conformations. By displacing protein structures along low-frequency normal modes and introducing mutations, MELT can either lock proteins into specific conformations or enhance dynamic behaviors along chosen normal modes.

## Key Features:
- Utilizes Normal Mode Analysis (NMA) to explore protein flexibility and dynamics.
- Introduces mutations that trap or favor dynamics, altering the protein's energy landscape.

## Installation
Prerequisites:

R for executing the normal mode analysis (NMA).
PyRosetta for structural relaxation and mutagenesis design.
Make sure to have the required libraries installed:

```
# R installation
install.packages("bio3d")
```

```
# PyRosetta (follow PyRosetta installation instructions)
pip install pyrosetta
```

## Step 1: Normal Mode Analysis (NMA)

## Parameters:
pdb: Input PDB structure (e.g., 2LZT).
mode: Normal mode used for displacement (e.g., mode 7).
direction: Displacement direction ("plus" or "minus").
scale: Scale of displacement (default is 10).

## Outputs:
Displaced PDB files: desloc_<mode>_<direction>.pdb

The displacement will generate altered conformations of the protein based on the selected normal mode, such as:

desloc_7_minus.pdb
desloc_7_plus.pdb

## Step 2: In Silico Mutagenesis and Relaxation
Once the protein structure is displaced, you can apply in silico mutagenesis using PyRosetta by running the Step2_DesignPyrosetta.py script.

Example Usage:

```
python Step2_DesignPyrosetta.py
```

Relaxation: After displacement, the structure is relaxed using the FastRelax protocol to minimize energy.

Mutagenesis: Selected regions of the protein are mutated to introduce the desired dynamics while maintaining structural integrity.

## Parameters:
pose: Displaced protein structure loaded into PyRosetta.
scorefxn: PyRosetta score function for energy calculations (e.g., ref2015_cart.wts).
chain: Chain of the protein to design (e.g., "A").
Outputs:
Relaxed and mutated PDB files:

desloc_7_relaxed.pdb

desloc_7_relaxed_designed.pdb

The mutated designs aim to either lock the protein into the displaced conformation or enhance the dynamics observed along the chosen normal mode.

## Workflow Summary
## Normal Mode Analysis:

Displace protein structures along low-frequency normal modes using Step1_NormalModes.R.

Output: Displaced PDB files.

In Silico Mutagenesis and Relaxation:

Use PyRosetta to apply mutagenesis and relax the displaced structures.

Output: Mutated, relaxed PDB files with desired dynamic properties.

Validate the protein dynamics using molecular dynamics, and if possible confirm the structures using AlphaFold
