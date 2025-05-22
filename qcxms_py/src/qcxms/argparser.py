import argparse
import sys
from dataclasses import fields

try:
    from .data import RunTypeData, WP
    from . import iomod # For check_progs, citation
    # from .utility import check_progs, citation # If they are moved there
except ImportError:
    # Fallbacks for standalone development
    from data import RunTypeData, WP
    class IomodMock: # Mock iomod for standalone testing
        @staticmethod
        def check_prog(pname, verbose=False, critical=True): print(f"Mock check_prog: {pname}")
    iomod = IomodMock()


def parse_arguments(argv_list=None) -> RunTypeData:
    """Parses command-line arguments and populates a RunTypeData object."""

    parser = argparse.ArgumentParser(
        description="QCxMS2: Calculate Electron Ionization Mass Spectra.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False # Custom help action will be added
    )

    # Create an instance of RunTypeData to get default values
    # These defaults will be used by argparse
    defaults = RunTypeData()

    # --- Special Help Arguments ---
    parser.add_argument(
        '-h', '--help', action='help', default=argparse.SUPPRESS,
        help='Show this help message and exit.'
    )
    # TODO: Implement advanced help later if needed, for now, all options are shown.
    # parser.add_argument(
    #     '--advanced', action='store_true', default=False,
    #     help='Show advanced options and exit.' # This would need custom action
    # )


    # --- Positional Argument ---
    parser.add_argument(
        'infile', type=str, 
        help='Input coordinate file in XYZ format (Angstroms).'
    )

    # --- Run Type Options ---
    runtype_group = parser.add_argument_group('Run Type Options')
    runtype_group.add_argument(
        '-ei', dest='mode_ei', action='store_true', default=(defaults.mode == 'ei'), # Default if mode is 'ei'
        help="Perform prediction of EI spectrum (default if no mode specified)."
    )
    runtype_group.add_argument(
        '-cid', dest='mode_cid', action='store_true', default=(defaults.mode == 'cid'),
        help="Perform prediction of CID spectrum."
    )
    # Note: Fortran code sets env.mode. Here we'll set it based on these flags later.
    # If both are true or false, need a priority or default. Default is 'ei'.

    runtype_group.add_argument(
        '-esim', dest='esim', action='store_true', default=defaults.esim,
        help="Simulate only different IEE distributions for a given reaction network from a previous calculation."
    )
    runtype_group.add_argument(
        '-oplot', dest='oplot', action='store_true', default=defaults.oplot,
        help="Plot only peaks from a previous QCxMS2 calculation (requires exp.dat as CSV)."
    )
    runtype_group.add_argument(
        '-ircrun', dest='ircrun', action='store_true', default=defaults.ircrun,
        help="Search in hessian calculation directory for an IRC and perform analysis."
    )
    runtype_group.add_argument(
        '-picktsrun', dest='picktsrun', action='store_true', default=defaults.picktsrun,
        help="Search reaction path for maxima/potential TSs from a previous calculation."
    )


    # --- General and Technical Options ---
    general_group = parser.add_argument_group('General and Technical Options')
    general_group.add_argument(
        '-chrg', type=int, default=defaults.chrg, dest='chrg',
        help="Charge of the input molecule."
    )
    general_group.add_argument(
        '-edist', type=str, default=defaults.edist, choices=['poisson', 'gaussian'], dest='edist',
        help="Distribution for Internal Electron Energy (IEE)."
    )
    general_group.add_argument(
        '-eimp0', type=WP, default=defaults.eimp0, dest='eimp0',
        help="Impact energy of electrons in eV for IEE distribution."
    )
    general_group.add_argument(
        '-eimpw', type=WP, default=defaults.eimpw, dest='eimpw',
        help="Width of IEE distribution."
    )
    general_group.add_argument(
        '-ieeatm', type=WP, default=defaults.ieeatm, dest='ieeatm',
        help="Energy per atom in eV for IEE distribution scaling."
    )
    general_group.add_argument(
        '-nots', dest='nots', action='store_true', default=defaults.nots,
        help="Take reaction energy (Delta E) instead of activation energy (Ea) -> no path search (quickmode for fragments)."
    )
    general_group.add_argument(
        '-T', type=int, default=defaults.cores, dest='cores',
        help="Total number of CPU cores to use for parallel calculations."
    )
    general_group.add_argument(
        '-nfrag', type=int, default=defaults.nfrag, dest='nfrag',
        help="Number of subsequent fragmentations to simulate."
    )
    general_group.add_argument(
        '-pthr', type=WP, default=defaults.pthr, dest='pthr',
        help="Intensity threshold (in %%) at which a fragment is selected for further fragmentation."
    )
    general_group.add_argument(
        '-tf', type=WP, default=defaults.tf, dest='tf',
        help="Time of flight in mass spectrometer, in microseconds (µs)."
    )
    general_group.add_argument(
        '-nsamples', type=int, default=defaults.nsamples, dest='nsamples',
        help="Number of simulated runs in Monte Carlo for intensity calculations."
    )
    general_group.add_argument(
        '-iseed', type=int, default=42, # Fortran default for iseed(1)
        help="Seed for random number generator. Default (42) or specific value."
    )


    # --- QC Calculation Options ---
    qc_group = parser.add_argument_group('Quantum Chemistry Calculation Options')
    qc_group.add_argument(
        '-geolevel', type=str, default=defaults.geolevel, dest='geolevel',
        help="Method for geometry optimization and path search (e.g., gfn2, gfn1, r2scan3c)."
    )
    qc_group.add_argument(
        '-tslevel', type=str, default=defaults.tslevel, dest='tslevel',
        help="Method for computing reaction barriers (e.g., gfn2, gfn1, r2scan3c)."
    )
    qc_group.add_argument(
        '-iplevel', type=str, default=defaults.iplevel, dest='iplevel',
        help="Method for computing IPs for charge assignment prescreening."
    )
    qc_group.add_argument(
        '-ip2level', type=str, default=defaults.ip2level if defaults.ip2level else None, # Handle empty string default
        help="Method for computing IPs for charge assignment refinement of critical cases (optional)."
    )
    qc_group.add_argument(
        '-qmprog', type=str, default=defaults.qmprog, choices=['orca', 'tm'], dest='qmprog',
        help="External QM program to use (orca or turbomole)."
    )
    qc_group.add_argument(
        '-tsfinder', type=str, default=defaults.tsfinder, choices=['neb', 'gsm', 'xtb'], dest='tsfinder',
        help="Pathfinder method (neb, gsm, or xtb internal)."
    )
    qc_group.add_argument(
        '-nebnormal', dest='nebnormal', action='store_true', default=defaults.nebnormal,
        help="Use normal (tighter) settings instead of loose settings for NEB search."
    )
    qc_group.add_argument(
        '-checkmult', dest='checkmult', action='store_true', default=defaults.checkmult,
        help="Check multiplicity of transition states."
    )
    qc_group.add_argument(
        '-pathmult', dest='pathmult', action='store_true', default=defaults.pathmult,
        help="Compute reaction path with the sum of multiplicities of products."
    )
    qc_group.add_argument(
        '-tsoptgmf', dest='tsoptgmf', action='store_true', default=defaults.tsoptgmf,
        help="Use special GMF optimizer in ORCA for TS optimization (requires ORCA dev version)."
    )
    qc_group.add_argument(
        '-tsoptact', dest='tsoptact', action='store_true', default=defaults.tsoptact,
        help="Optimize TS in ORCA by specifying active atoms."
    )
    qc_group.add_argument(
        '-notsgeo', dest='tsgeodesic_false', action='store_true', default=not defaults.tsgeodesic,
        help="Do NOT use geodesic interpolation as guess for NEB/GSM runs."
    )
    qc_group.add_argument(
        '-tsnodes', type=int, default=defaults.tsnds, dest='tsnds',
        help="Number of nodes for path search (NEB/GSM)."
    )
    qc_group.add_argument(
        '-fermi', dest='fermi', action='store_true', default=defaults.fermi,
        help="Apply Fermi smearing in SCF calculations."
    )
    qc_group.add_argument(
        '-exstates', type=int, default=defaults.exstates, dest='exstates',
        help="Number of excited states to include via TDDFT (requires DFT functional)."
    )


    # --- CREST msreact Options ---
    msreact_group = parser.add_argument_group('CREST msreact Options')
    msreact_group.add_argument(
        '-msnoiso', dest='msnoiso', action='store_true', default=defaults.msnoiso,
        help="Sort out non-dissociated (isomeric) fragments in CREST msreact."
    )
    msreact_group.add_argument(
        '-msiso', dest='msiso', action='store_true', default=defaults.msiso,
        help="Sort out dissociated fragments in CREST msreact (i.e., only consider isomers)."
    )
    # msfulliso is default True in Fortran. Argparse action='store_false' if flag means disable.
    msreact_group.add_argument(
        '--no-msfulliso', dest='msfulliso', action='store_false', default=defaults.msfulliso, 
        help="Disable rearrangements for subsequent fragmentations (default is to allow them)."
    )
    msreact_group.add_argument(
        '-msnbonds', type=int, default=defaults.msnbonds, dest='msnbonds',
        help="Max number of bonds between atom pairs for applying repulsive potential in CREST msreact."
    )
    msreact_group.add_argument(
        '-msnshifts', type=int, default=defaults.msnshifts, dest='msnshifts',
        help="Perform N optimizations with randomly shifted atom positions in CREST msreact."
    )
    msreact_group.add_argument(
        '-msnshifts2', type=int, default=defaults.msnshifts2, dest='msnshifts2',
        help="Perform N optimizations with randomly shifted atom positions and repulsive potential in CREST msreact."
    )
    # msmolbar is default True.
    msreact_group.add_argument(
        '--no-msmolbar', dest='msmolbar', action='store_false', default=defaults.msmolbar,
        help="Disable sorting out duplicates by molbar codes (requires molbar)."
    )
    msreact_group.add_argument(
        '-msinchi', dest='msinchi', action='store_true', default=defaults.msinchi,
        help="Sort out duplicates by InChI codes (requires obabel)."
    )
    msreact_group.add_argument(
        '-msfragdist', type=WP, default=defaults.msfragdist, dest='msfragdist',
        help="Separate fragments from each other before TS search (distance in Angstrom) in CREST msreact."
    )
    msreact_group.add_argument(
        '--mskeepdir', dest='mskeepdir', action='store_true', default=defaults.mskeepdir, # Fortran default was .false.
        help="Keep the MSDIR directory with CREST msreact constrained optimizations."
    )

    # --- Thermochemistry and Rate Calculation Options ---
    thermo_group = parser.add_argument_group('Thermochemistry and Rate Calculation Options')
    thermo_group.add_argument(
        '--nobhess', dest='bhess', action='store_false', default=defaults.bhess,
        help="Do NOT add thermochemical contribution (ZPE/thermal) to activation energy."
    )
    thermo_group.add_argument(
        '--usetemp', dest='notemp', action='store_false', default=defaults.notemp,
        help="Use G(RRHO) thermal corrections instead of only ZPVE (default is ZPVE only)."
    )
    thermo_group.add_argument(
        '-atmassh', type=WP, default=defaults.hmass, dest='hmass',
        help="Modify mass of Hydrogen atoms (e.g., for scaling R-H frequencies)."
    )
    thermo_group.add_argument(
        '-scaleeinthdiss', type=WP, default=defaults.scaleeinthdiss, dest='scaleeinthdiss',
        help="Scaling factor for internal energy in -H or -H2 abstractions."
    )
    thermo_group.add_argument(
        '-scaleker', type=WP, default=defaults.scaleker, dest='scaleker',
        help="Scaling factor for Kinetic Energy Release (KER) upon fragmentation."
    )
    thermo_group.add_argument(
        '-tfscale', type=WP, default=defaults.tfscale, dest='tfscale',
        help="Scaling factor for time-of-flight of subsequent fragmentations."
    )
    thermo_group.add_argument(
        '--noeatomscale', dest='eatomscale', action='store_false', default=defaults.eatomscale,
        help="Do NOT scale IEE of subsequent fragmentations by atom number ratio."
    )
    thermo_group.add_argument(
        '-sumreacscale', type=WP, default=defaults.sumreacscale, dest='sumreacscale',
        help="Scaling factor for sum_reac in Monte Carlo simulations."
    )
    thermo_group.add_argument(
        '-erelaxtime', type=WP, default=defaults.erelaxtime * 1e12, # Convert ps to input unit for user
        help="Assumed relaxation time in picoseconds for H-dissociation scaling."
    )
    thermo_group.add_argument(
        '--noKER', dest='calcKER', action='store_false', default=defaults.calcKER,
        help="Do NOT compute Kinetic Energy Release (KER)."
    )
    thermo_group.add_argument(
        '--noreoptts', dest='reoptts', action='store_false', default=defaults.reoptts,
        help="Do NOT reoptimize Transition States after path search."
    )
    thermo_group.add_argument(
        '-hotip', dest='hotip', action='store_true', default=defaults.hotip,
        help="Use unoptimized (hot) fragment geometries for IP calculations."
    )
    thermo_group.add_argument(
        '-sortoutcascade', dest='sortoutcascade', action='store_true', default=defaults.sortoutcascade,
        help="Sort out reaction pathways with more than one maximum (cascade reactions)."
    )
    thermo_group.add_argument(
        '-rrkm', dest='eyring_false', action='store_true', default=not defaults.eyring, # if true, eyring becomes false
        help="Use simplified RRKM equation instead of Eyring for rate constants."
    )
    thermo_group.add_argument(
        '-eyzpve', dest='eyzpve', action='store_true', default=defaults.eyzpve,
        help="Use Eyring equation with ZPVE only (Delta H like) instead of full Delta G for rates."
    )
    thermo_group.add_argument(
        '-sthr', type=int, default=defaults.sthr, dest='sthr',
        help="RRHO vibrational frequency cutoff (cm-1) for entropy calculations."
    )
    thermo_group.add_argument(
        '-ithr', type=int, default=defaults.ithr, dest='ithr',
        help="Imaginary frequency cutoff (cm-1) for RRHO vibrational calculations."
    )
    thermo_group.add_argument(
        '-nthermosteps', type=int, default=defaults.nthermosteps, dest='nthermosteps',
        help="Number of increments to compute thermal corrections for IEE distribution."
    )

    # --- Plotting Options ---
    plot_group = parser.add_argument_group('Plotting Options')
    plot_group.add_argument(
        '-plotk', dest='plot', action='store_true', default=defaults.plot,
        help="Plot k(E) curves for analysis."
    )
    plot_group.add_argument(
        '--noisotope', dest='noiso', action='store_true', default=defaults.noiso,
        help="Do NOT plot isotope patterns in the mass spectrum."
    )
    plot_group.add_argument(
        '-int_masses', dest='int_masses', action='store_true', default=defaults.int_masses,
        help="Plot m/z values as integers only."
    )
    plot_group.add_argument(
        '-mthr', type=WP, default=defaults.mthr, dest='mthr',
        help="Minimum m/z threshold for plotting peaks."
    )
    plot_group.add_argument(
        '-intthr', type=WP, default=defaults.intthr, dest='intthr',
        help="Minimum peak intensity threshold (in %%) for plotting."
    )
    plot_group.add_argument(
        '-noplotisos', dest='noplotisos', action='store_true', default=defaults.noplotisos,
        help="Do NOT plot isomers in the mass spectrum (assumes they fragment further)."
    )


    # --- Special/Debugging Options ---
    special_group = parser.add_argument_group('Special and Debugging Options')
    special_group.add_argument(
        '-dxtbparam', dest='dxtb', action='store_true', default=defaults.dxtb,
        help="Use special dxtb parameters (requires dxtb_param.txt in starting directory)."
    )
    special_group.add_argument(
        '-noirc', dest='noirc', action='store_true', default=defaults.noirc,
        help="Set all imaginary frequencies to 100 cm-1 (effective in -esim mode only)."
    )
    special_group.add_argument(
        '-cneintscale', dest='cneintscale', action='store_true', default=defaults.cneintscale,
        help="Scale internal energy of subsequent fragmentations by coordination number of reaction atoms."
    )
    special_group.add_argument(
        '-printlevel', type=int, default=defaults.printlevel, choices=[1, 2, 3], dest='printlevel',
        help="Printout verbosity level (1: normal, 2: verbose, 3: debug)."
    )
    special_group.add_argument(
        '-logo', dest='logo', action='store_true', default=defaults.logo,
        help="Print QCxMS2 logo at startup."
    )
    special_group.add_argument(
        '--keepdir', dest='removedirs_false', action='store_true', default=not defaults.removedirs,
        help="Do NOT remove temporary fragment directories after calculation (for diagnostics)."
    )
    special_group.add_argument(
        '-ignoreip', dest='ignoreip', action='store_true', default=defaults.ignoreip,
        help="Ignore IPs; both fragments get full intensity (for testing)."
    )
    special_group.add_argument(
        '-chargeatomscale', dest='chargeatomscale', action='store_true', default=defaults.chargeatomscale,
        help="Larger fragment gets the charge, IPs are ignored (for testing)."
    )
    special_group.add_argument(
        '-fixe', type=WP, default=defaults.fixe, dest='fixe',
        help="Simulate only one specific energy (eV) for IEE (for testing; 0 means use distribution)."
    )
    special_group.add_argument(
        '-topocheck', type=str, default=defaults.topocheck, dest='topocheck',
        help="Method for topology check after optimization ('molbar', 'inchi', or empty string to disable)."
    )

    # --- CID Specific Options ---
    cid_group = parser.add_argument_group('CID Specific Options (Highly Experimental)')
    cid_group.add_argument(
        '-protonate', dest='prot', action='store_true', default=defaults.prot,
        help="Perform search for favored protomers (for ESI-MS like CID setup)."
    )
    cid_group.add_argument(
        '-deprotonate', dest='deprot', action='store_true', default=defaults.deprot,
        help="Perform search for favored deprotonated structures (for ESI-MS like CID setup)."
    )
    # CID mode is set by -cidauto, -cidtemprun, -cidforced
    cid_group.add_argument(
        '-cidauto', dest='cid_mode_auto', action='store_true', default=(defaults.cid_mode == 1 and defaults.mode == 'cid'),
        help="Auto mode for CID simulation (default if -cid is active)."
    )
    cid_group.add_argument(
        '-cidtemprun', dest='cid_mode_temprun', action='store_true', default=(defaults.cid_mode == 2 and defaults.mode == 'cid'),
        help="CID simulation without collisions (only thermal heating in ESI simulation)."
    )
    cid_group.add_argument(
        '-cidforced', dest='cid_mode_forced', action='store_true', default=(defaults.cid_mode == 3 and defaults.mode == 'cid'),
        help="CID simulation with only collisions (no thermal heating in ESI simulation)."
    )
    cid_group.add_argument(
        '-elab', type=WP, default=defaults.cid_elab, dest='cid_elab',
        help="Collision energy in laboratory frame (eV) for CID."
    )
    cid_group.add_argument(
        '-esi', type=WP, default=defaults.cid_esi, dest='cid_esi',
        help="Internal energy scaling in ESI (eV)."
    )
    cid_group.add_argument(
        '-esiatom', type=WP, default=defaults.cid_esiatom, dest='cid_esiatom',
        help="Internal energy scaling per atom in ESI (eV/atom)."
    )
    cid_group.add_argument(
        '-esiw', type=WP, default=defaults.cid_esiw, dest='cid_esiw',
        help="Width of internal energy distribution in ESI (eV)."
    )
    cid_group.add_argument(
        '-collw', type=WP, default=defaults.cid_collw, dest='cid_collw',
        help="Width of collision energy distribution (eV)."
    )
    cid_group.add_argument(
        '-maxcoll', type=int, default=defaults.cid_maxcoll, dest='cid_maxcoll',
        help="Maximum number of collisions in CID."
    )
    cid_group.add_argument(
        '-lchamb', type=WP, default=defaults.cid_lchamb, dest='cid_lchamb',
        help="Length of collision chamber (m)."
    )
    cid_group.add_argument(
        '-mgas', type=WP, default=defaults.cid_mgas, dest='cid_mgas',
        help="Mass of collision gas (amu), e.g., Argon = 39.948."
    )
    cid_group.add_argument(
        '-rgas', type=WP, default=defaults.cid_rgas, dest='cid_rgas',
        help="Van der Waals radius of collision gas (Bohr)."
    )
    cid_group.add_argument(
        '-pgas', type=WP, default=defaults.cid_pgas, dest='cid_pgas',
        help="Pressure of collision gas (Pa)."
    )
    cid_group.add_argument(
        '-cidscool', type=WP, default=defaults.cid_scool, dest='cid_scool',
        help="Scaling factor for collisional cooling of ions in CID mode."
    )
    cid_group.add_argument(
        '-solv', dest='solv', action='store_true', default=defaults.solv,
        help="Use solvation model for barrier and energy calculations."
    )
    
    # --- Parse Arguments ---
    if argv_list is None: # Allow testing with a list of strings
        argv_list = sys.argv[1:]
    
    # Handle empty argv_list: if called with no args, print help and exit
    if not argv_list:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(argv_list)

    # --- Populate RunTypeData ---
    env = RunTypeData()

    # Transfer parsed values to RunTypeData instance
    # For most, dest matches attribute name, so vars(args) can be used.
    for f in fields(RunTypeData):
        if hasattr(args, f.name):
            setattr(env, f.name, getattr(args, f.name))
    
    # Handle special cases / logic from Fortran
    env.infile = args.infile # Positional argument

    # Mode selection (EI/CID)
    if args.mode_cid or args.cid_mode_auto or args.cid_mode_temprun or args.cid_mode_forced:
        env.mode = 'cid'
        if args.cid_mode_temprun:
            env.cid_mode = 2
            env.cid_maxcoll = 0 # As per Fortran logic
        elif args.cid_mode_forced:
            env.cid_mode = 3
        else: # auto or just -cid
            env.cid_mode = 1
    elif args.mode_ei: # Explicit -ei or by default if no CID flags
        env.mode = 'ei'
    # If neither mode_cid nor mode_ei is set by user, it defaults to 'ei' from RunTypeData init.

    if args.tsgeodesic_false: # Flag means set to false
        env.tsgeodesic = False
    
    if args.eyring_false: # Flag means set eyring to false (i.e. use RRKM)
        env.eyring = False
        
    if args.removedirs_false: # Flag means set removedirs to false
        env.removedirs = False
    
    # Convert erelaxtime from ps back to seconds
    env.erelaxtime = env.erelaxtime * 1.0e-12

    # Seed for random number generator (currently not used directly in RunTypeData for generation)
    # The Fortran code had `iseed(1)` for this. We'll store it if needed by other parts.
    # For now, `args.iseed` is just stored if it's added to RunTypeData.
    # If `RunTypeData` needs an `iseed` field:
    # env.iseed = [args.iseed] # or handle as a list if that's the plan

    return env


def check_settings_py(env: RunTypeData):
    """Prints a summary of important settings. Python equivalent of Fortran check_settings."""
    print("Settings:")
    print("*" * 60)
    if env.chrg == 1: print("           + Positive Ion mode +")
    elif env.chrg == -1: print("           - Negative Ion mode -")
    elif env.chrg > 1 or env.chrg < -1:
        print(f"Charge of input molecule was set to: {env.chrg}")
        print("!!!Warning multiply charged ions not yet tested!!!")

    print(f"IEE distribution: {env.edist}")
    print(f"Width of IEE distribution (eimpw) is: {env.eimpw:.2f}")
    print(f"Energy per atom (ieeatm) is: {env.ieeatm:.2f} eV")
    print(f"eimp0 is: {env.eimp0:.2f} eV")

    if env.exstates > 0 and env.tslevel not in ['wb97x3c', 'r2scan3c', 'pbeh3c']: # Add other DFT functionals if supported
        print("Warning: Excited states via TD-DFT only possible for DFT functionals")
        # sys.exit("Exiting due to incompatible TD-DFT settings.") # Or handle as error

    if env.mode == "ei":
        print("Spectral mode: EI-MS")
    elif env.mode == "cid":
        if env.cid_mode == 1: print("Spectral mode: CID auto (!WARNING: highly experimental!)")
        elif env.cid_mode == 2: print("Spectral mode: CID temprun (no collisions simulated)")
        elif env.cid_mode == 3:
            print("Spectral mode: CID forced collisions (!WARNING: not yet implemented)")
            print("Aborting...")
            sys.exit(1)
    
    if env.chrg < 0: print("Negative ion mode chosen, Hybrid DFT with diffuse basis sets recommended!")
    if env.qmprog != "orca": print("Warning! Currently only ORCA is fully supported as external QM program.")

    print(f"Path search method: {env.tsfinder}")
    if env.tsfinder == "gsm":
        print("Using double-ended growing string methods for path search.")
        if env.tsnds == 9: env.tsnds = 15 # Default adjustment from Fortran
        print(f"Number of nodes for GSM set to: {env.tsnds}")
    
    if env.geolevel == "gxtb" and env.tsfinder == 'neb':
        print("Warning: Special xTB version for NEB search might be necessary.")

    print("QM methods used:")
    print(f"  Level for geometry optimizations and path searches: {env.geolevel}")
    print(f"  Level for reaction energy and barrier calculations: {env.tslevel}")
    print(f"  Level for IP prescreening: {env.iplevel}")
    if env.ip2level and env.ip2level != env.iplevel:
        print(f"  Level for IP refinement: {env.ip2level}")

    print("Runtime settings used:")
    print(f"  Number of cores used: {env.cores}")
    print(f"  Number of allowed subsequent fragmentations: {env.nfrag}")
    print(f"  Intensity threshold for further fragmentation: {env.pthr:.2f}%")
    if env.exstates > 0: print(f"  Number of excited states to include via TDDFT: {env.exstates}")
    if env.solv: print("Solvation selected for barrier and energy calculations.")
    print("*" * 60)


def check_programs_py(env: RunTypeData):
    """Checks for necessary external programs. Python equivalent of Fortran check_progs."""
    print("External programs used:")
    iomod.check_prog('xtb', verbose=True, critical=True)
    iomod.check_prog('crest', verbose=True, critical=True)

    if env.reoptts or env.tslevel not in ["gfn2", "gfn1", ""] or \
       env.iplevel not in ["gfn2", "gfn1", ""] or \
       (env.ip2level and env.ip2level not in ["gfn2", "gfn1", ""]):
        iomod.check_prog('orca', verbose=True, critical=True)
    
    if env.tsfinder == "gsm":
        iomod.check_prog('gsm.orca', verbose=True, critical=True)
        iomod.check_prog('tm2orca.py', verbose=True, critical=True)

    xtb_path = shutil.which('xtb')
    if xtb_path:
        os.environ['XTBEXE'] = xtb_path
        print(f"XTB path set for ORCA: {xtb_path}")
    else: # Already checked by check_prog, but good to be explicit if XTBEXE is crucial
        print("Warning: XTB not found, XTBEXE environment variable not set.")

    if env.msmolbar or env.topocheck == "molbar":
        iomod.check_prog('molbar', verbose=True, critical=True)
    if env.msinchi or env.topocheck == "inchi":
        iomod.check_prog('obabel', verbose=True, critical=True) # InChI usually via obabel

    if env.tsgeodesic:
        iomod.check_prog('geodesic_interpolate', verbose=True, critical=False) # Not always critical

    levels_to_check_gxtb = [env.tslevel, env.iplevel, env.ip2level, env.geolevel]
    if any(level == 'gxtb' for level in levels_to_check_gxtb):
        iomod.check_prog('gxtb', verbose=True, critical=True)
        # Check for gxtb specific files - adapt path as needed
        # home_dir = Path.home()
        # if not (home_dir / ".gxtb").exists(): sys.exit("Error: ~/.gxtb file required for gxtb not found.")
        # if not (home_dir / ".ceh").exists(): sys.exit("Error: ~/.ceh file required for gxtb not found.")
        # if not (home_dir / ".basisq").exists(): sys.exit("Error: ~/.basisq file required for gxtb not found.")
        print("Note: gxtb usage implies ~/.gxtb, ~/.ceh, ~/.basisq files are required.")

    print("*" * 60 + "\n")


if __name__ == '__main__':
    # Example of how to use the parser
    # To test, run: python argparser.py test.xyz -T 8 --no-msfulliso -eimp0 65.0
    # Or: python argparser.py test.xyz (to see defaults applied)
    
    if len(sys.argv) == 1: # If script is called without any arguments
        # Create a dummy test.xyz file for the default run
        Path("test.xyz").write_text("1\nH 0 0 0\n")
        print("Running with default test file 'test.xyz' and minimal options:")
        test_argv = ["test.xyz"]
        env_settings = parse_arguments(test_argv)
    else:
        env_settings = parse_arguments() # Uses sys.argv[1:]

    print("\n--- Parsed RunTypeData ---")
    import json
    # A simple way to see all settings, though dataclasses.asdict would be better if not for WP type
    # For now, just print a few key ones or iterate through fields
    # print(json.dumps(env_settings.__dict__, indent=2))
    for f in fields(RunTypeData):
        print(f"{f.name}: {getattr(env_settings, f.name)}")
    
    print("\n--- Checking Settings & Programs ---")
    check_settings_py(env_settings)
    check_programs_py(env_settings)
    # citation_py() # If citation is moved here

    if Path("test.xyz").exists() and len(sys.argv) == 1:
        Path("test.xyz").unlink()
```
