library(bio3d)

# Read in a PDB file for a protein structure
pdb <- read.pdb("2LZT")

# Trim the PDB file to remove certain atoms:
# Exclude HETATM entries and hydrogen atoms
# The 'operator = "OR"' indicates a logical OR operation for selection
# 'inverse = TRUE' means all atoms except those matching the criteria will be kept
pdb <- trim.pdb(pdb, atom.select(pdb, type = 'HETATM', elesy = "H", 
                                 operator = "OR", inverse = TRUE))

# Perform normal mode analysis (NMA) on the trimmed structure
# 'outmodes = "noh"' indicates excluding hydrogen atoms from output modes
# 'rm.wat = TRUE' removes water molecules from the analysis
nma <- aanma(pdb, outmodes = 'noh', rm.wat = TRUE)

# Function to apply displacement in a specified direction and mode
# Args:
#   pdb: The original PDB object
#   nma: The NMA results containing the mode vectors
#   mode: The mode number to use for displacement
#   direction: The direction of displacement, either "plus" or "minus"
#   scale: The scaling factor for the displacement (default is 10)
#   output_file: The name of the output PDB file
apply_displacement <- function(pdb, nma, mode, direction, scale = 10) {
  # Copy the original pdb object to preserve it
  pdb_modified <- pdb
  
  # Apply displacement based on the specified direction
  if (direction == "minus") {
    pdb_modified$xyz <- pdb_modified$xyz - scale * nma$modes[, mode]
  } else if (direction == "plus") {
    pdb_modified$xyz <- pdb_modified$xyz + scale * nma$modes[, mode]
  } else {
    stop("Invalid direction. Please specify 'plus' or 'minus'.")
  }
  
  # Write the modified PDB structure to a file
  write.pdb(pdb_modified, file = paste("desloc_",mode,"_",direction,".pdb", sep = ""))
}

# Example usage: Displacing using mode 7 in the minus direction
apply_displacement(pdb, nma, mode = 7, direction = "minus", 
                   scale = 10)

# Example usage: Displacing using mode 7 in the plus direction
apply_displacement(pdb, nma, mode = 7, direction = "plus", 
                   scale = 10)