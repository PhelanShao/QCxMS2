# QCxMS_Py: Python Port of QCxMS2

QCxMS_Py is a Python port of the QCxMS2 program, which is designed for the calculation and simulation of Electron Ionization (EI) mass spectra and Collision-Induced Dissociation (CID) mass spectra. The program automates reaction network discovery and employs quantum chemistry calculations to predict mass spectra.

This Python version aims to replicate and extend the functionality of the original Fortran-based QCxMS2, leveraging Python's extensive scientific computing ecosystem.

## Features (Planned/In Progress)

*   Generation of fragmentation pathways.
*   Quantum mechanical calculations for energies, geometries, and properties of molecules and fragments (interfacing with XTB, ORCA).
*   Calculation of Ionization Potentials (IPs) and statistical charge assignment.
*   Simulation of isotopic patterns.
*   Simulation of Collision-Induced Dissociation (CID) processes.
*   Monte Carlo simulation of fragmentation kinetics to produce mass spectra.
*   Plotting of final mass spectra.

## Project Structure

*   `src/qcxms/`: Main Python source code for the QCxMS library.
    *   `main.py`: Main executable script.
    *   `argparser.py`: Command-line argument parsing.
    *   `data.py`: Core data structures (e.g., `RunTypeData`).
    *   `iomod.py`: File I/O and system interaction utilities.
    *   `utility.py`: General utility functions.
    *   `qmmod.py`: Interface for QM calculations (XTB, ORCA).
    *   `charges.py`: Charge assignment logic.
    *   `isotopes.py`: Isotopic pattern generation.
    *   `iee.py`: Impact Excess Energy (IEE) distribution generation.
    *   `fragmentation.py`: Fragmentation pathway generation and orchestration.
    *   `reaction.py`: Reaction energy calculation and product processing.
    *   `cid.py`: CID simulation logic.
    *   `structools.py`: Molecular structure analysis tools.
    *   `rmsd.py`: RMSD calculation.
    *   `boxmuller.py`: Box-Muller transform for random number generation.
    *   *(Placeholder for `tsmod.py`, `mcsimu.py`, `plot.py`)*
*   `tests/`: Unit and integration tests (to be developed).
*   `pyproject.toml`: Project metadata, dependencies, and build configuration.
*   `README.md`: This file.
*   `.gitignore`: Specifies intentionally untracked files that Git should ignore.
*   `.github/workflows/`: GitHub Actions CI/CD workflows.

## Installation

This project uses a `pyproject.toml` file and can be built and installed using standard Python packaging tools like `pip` and `build`.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-repo/qcxms_py.git # Replace with actual URL
    cd qcxms_py
    ```

2.  **Create a virtual environment (recommended):**
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # On Windows use `.venv\Scripts\activate`
    ```

3.  **Install the package:**
    For development (editable install):
    ```bash
    pip install -e .
    ```
    For a regular install:
    ```bash
    pip install .
    ```
    This will also install dependencies like NumPy and TOML.

## Usage

The program can be run using the installed script entry point:

```bash
qcxms_run input.xyz [options]
```

Replace `input.xyz` with your input molecular geometry file and `[options]` with any desired command-line flags.

To see available options:
```bash
qcxms_run --help
```

## Dependencies

Core dependencies are listed in `pyproject.toml` and include:
*   Python (>=3.8)
*   `numpy`
*   `toml`

External QM software is also required for full functionality:
*   XTB (for semi-empirical calculations)
*   ORCA (for DFT and other QM calculations)
*   CREST (for fragment generation using `--msreact`)
*   Optional: `molbar` or `obabel` for certain topology checks.

Ensure these programs are installed and available in your system's PATH.

## Contributing

Contributions are welcome! Please refer to the (forthcoming) `CONTRIBUTING.md` for guidelines.

## License

This project is licensed under the Apache-2.0 License. See the `pyproject.toml` file for details.
(The original QCxMS2 and its components like mctc-lib, PlotMS, etc., may have their own licenses - typically GPL or LGPL. This Python port is intended as a new implementation.)
```
