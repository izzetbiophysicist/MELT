library(bio3d)

# Read in a PDB file for a protein structure
pdb <- read.pdb("2LZT")

# Trim the PDB file to remove certain atoms:
# Exclude HETATM entries and hydrogen atoms
pdb <- trim.pdb(pdb, atom.select(pdb, type = 'HETATM', elesy = "H", 
                                 operator = "OR", inverse = TRUE))

# Perform normal mode analysis (NMA) on the trimmed structure
nma <- aanma(pdb, outmodes = 'noh', rm.wat = TRUE)

# Function to apply displacement in a specified manner (Factor or Geometrical)
# Args:
#   pdb: The original PDB object
#   nma: The NMA results containing the mode vectors
#   mode: The mode number to use for displacement
#   direction: The direction of displacement, either "plus" or "minus"
#   scale: The scaling factor for the Factor method (default is 10)
#   rms: The RMS value for the Geometrical method
#   selection: Atom selection for the Geometrical method ('calpha' or 'all')
#   method: The method to use for displacement ('Factor' or 'Geometrical')
apply_displacement <- function(pdb, nma, mode, direction, scale = 10, rms = 1.0, 
                               selection = 'calpha', method = 'Factor') {
  if (method == 'Factor') {
    # Displacement using the Factor method
    pdb_modified <- pdb

    # Apply displacement based on the specified direction
    if (direction == "minus") {
      pdb_modified$xyz <- pdb_modified$xyz - scale * nma$modes[, mode]
    } else if (direction == "plus") {
      pdb_modified$xyz <- pdb_modified$xyz + scale * nma$modes[, mode]
    } else {
      stop("Invalid direction. Please specify 'plus' or 'minus'.")
    }

    # Generate output file name
    output_file <- paste("desloc_", mode, "_", direction, "_factor.pdb", sep = "")

    # Write the modified PDB structure to a file
    write.pdb(pdb_modified, file = output_file)
  } else if (method == 'Geometrical') {
    # Displacement using the Geometrical method
    xyz <- pdb$xyz

    # Select atom indices based on input
    if (selection == 'calpha') {
      inds <- atom.select(pdb, elety = 'CA')$xyz
    } else if (selection == 'all') {
      inds <- atom.select(pdb, string = 'protein')$xyz
    } else {
      stop("Invalid selection. Choose 'calpha' or 'all'.")
    }

    # Compute scaling factor
    vec <- nma$modes[, mode]
    mag <- sqrt(length(inds)) * rms / sqrt(sum(vec^2))

    # Apply displacement
    displaced_xyz <- xyz
    displaced_xyz[inds] <- xyz[inds] + vec[inds] * mag

    # Generate output file name
    output_file <- paste("desloc_", mode, "_", direction, "_geometrical.pdb", sep = "")

    # Write to PDB
    write.pdb(pdb = pdb, xyz = displaced_xyz, file = output_file)
  } else {
    stop("Invalid method. Please specify 'Factor' or 'Geometrical'.")
  }
}

# Example usage: Displacing using mode 7 in the minus direction with Factor method
apply_displacement(pdb, nma, mode = 7, direction = "minus", 
                   scale = 10, method = 'Factor')

# Example usage: Displacing using mode 7 in the plus direction with Geometrical method
apply_displacement(pdb, nma, mode = 7, direction = "plus", 
                   rms = 1.0, selection = 'calpha', method = 'Geometrical')
