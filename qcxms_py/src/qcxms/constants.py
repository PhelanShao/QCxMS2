"""
Constants used throughout the QCxMS package.
All physical constants and conversion factors are defined here to avoid duplication.
"""

import math

# Precision parameter
WP = float

# Physical constants
BOHR = 0.52917726  # Bohr radius in Angstrom
AUTOEV = 27.211386245988  # Hartree to eV conversion
EVTOKCAL = 23.060547830618307  # eV to kcal/mol conversion
HARTREE_TO_EV = AUTOEV  # Alias for clarity
KB_EV_K = 8.6173332621415e-5  # Boltzmann constant in eV/K
PI = math.pi
PI_CONST = PI  # Alias for compatibility

# Conversion factors
EVTOAU = 1.0 / AUTOEV  # eV to Hartree
AUTOKCAL = EVTOKCAL  # Hartree to kcal/mol (via eV)

# Default values for calculations
DEFAULT_TEMP = 298.15  # Default temperature in Kelvin
DEFAULT_PRESSURE = 1.0  # Default pressure in atm

# Thresholds and cutoffs
DEFAULT_STHR = 150  # Default RRHO cutoff in cm-1
DEFAULT_ITHR = 100  # Default imaginary mode cutoff in cm-1
DEFAULT_PTHR = 1.0  # Default intensity threshold in %

# File extensions and names
XYZ_EXT = ".xyz"
DAT_EXT = ".dat"
OUT_EXT = ".out"

# Program names for external calls
PROGRAM_XTB = "xtb"
PROGRAM_CREST = "crest"
PROGRAM_ORCA = "orca"
PROGRAM_MOLBAR = "molbar"
PROGRAM_OBABEL = "obabel"