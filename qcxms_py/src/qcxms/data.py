from dataclasses import dataclass, field
import time
import sys

# Precision parameter
WP = float

# Bohr constant
BOHR: WP = 0.52917726

@dataclass
class Timer:
    times: int = 0
    rate: WP = 0.0  # System clock rate
    t: list[list[WP]] = field(default_factory=list)  # Stores timing data
    names: list[str] = field(default_factory=list)   # Stores names of timers

    def __post_init__(self):
        # Initialize t and names based on times if needed, or handle dynamically
        # For now, methods will handle dynamic allocation as in Fortran
        if self.times > 0 and not self.t:
            self.t = [[0.0] * 3 for _ in range(self.times)]
        if self.times > 0 and not self.names:
            self.names = [""] * self.times
        
        # Get the system clock rate (frequency)
        # In Python, time.perf_counter() returns seconds, so rate is effectively 1.0
        # However, to mimic the Fortran logic of system_clock(count_rate=rate),
        # we can use time.get_clock_info('perf_counter').resolution or simply 1.0
        # For simplicity, we'll assume perf_counter gives high enough resolution
        # and the rate calculation t(n,3) / rate is handled by direct differences.
        # The Fortran 'rate' is more about ticks per second.
        # Python's time.perf_counter() returns seconds, so differences are already in seconds.
        pass


    def init(self, n: int):
        self.times = n
        # self.rate is tricky. In Fortran, system_clock gives integer counts and a rate.
        # In Python, time.perf_counter() gives float seconds.
        # We don't strictly need 'rate' if we directly use perf_counter differences.
        self.t = [[0.0] * 3 for _ in range(n)] # [[start, stop, accumulated_duration], ...]
        self.names = [""] * n

    def clear(self):
        self.t = []
        self.names = []
        self.times = 0

    def start(self, n: int, inp: str):
        if 0 <= n < self.times:
            self.names[n] = inp
            self.t[n][0] = time.perf_counter()
        else:
            # Handle error or dynamic resizing if necessary
            # For now, assume n is always valid as per Fortran fixed size post-init
            print(f"Error: Timer index {n} out of bounds.")


    def stop(self, n: int):
        if 0 <= n < self.times:
            self.t[n][1] = time.perf_counter()
            self.t[n][2] += (self.t[n][1] - self.t[n][0])
        else:
            # Handle error
            print(f"Error: Timer index {n} out of bounds.")


@dataclass
class RunTypeData:
    # --- GENERAL data
    chrg: int = 1
    nfrags: int = 1
    gsmnodes: int = 5
    nsamples: int = 100000
    pthr: WP = 1.0
    temp: int = 298
    tsnds: int = 9
    threads: int = 4
    cores: int = 4

    # --- various names and flags
    infile: str = "" # Max length 80
    geolevel: str = 'gfn2' # Max length 80
    iplevel: str = 'gfn2' # Max length 80
    ip2level: str = 'gfn2' # Max length 80
    tslevel: str = 'gfn2' # Max length 80
    pathlevel: str = 'none' # Max length 80

    # external programs
    qmprog: str = 'orca' # Max length 80
    tsfinder: str = 'neb' # Max length 80
    topocheck: str = 'molbar' # Max length 80
    tsgeodesic: bool = True
    nebnormal: bool = False
    tsoptgmf: bool = False
    
    # some global parameters
    path: str = "" # Max length 1024, TODO FIXME can be removed
    index: str = "" # Max length 1024, TODO FIXME can be removed
    startdir: str = "" # Max length 1024
    
    # global run data
    # fraglist: list[str] = field(default_factory=lambda: [""] * 10000) # list of all fragment indices
    restartrun: bool = False

    # --- various parameter
    mode: str = 'ei' # Max length 10, "EI" or "CID"
    edist: str = 'poisson' # Max length 80, 'gaussian' or 'poisson'
    eimp0: WP = 70.0
    eimpw: WP = 0.1
    ieeatm: WP = 0.8
    fixe: WP = 0.0
    tf: WP = 50.0 # microseconds
    ehomo: WP = 0.0 # To be initialized or calculated
    hmass: WP = 1.0801
    scaleeinthdiss: WP = 0.5
    scaleker: WP = 1.0
    erelaxtime: WP = 5.0e-12
    sumreacscale: WP = 1.0
    eatomscale: bool = True
    chargeatomscale: bool = False
    exstates: int = 0

    # --- CREST MSREACT settings
    msnoiso: bool = False
    msfulliso: bool = True
    msiso: bool = False
    msnbonds: int = 3
    msnshifts: int = 200
    msnshifts2: int = 0
    msinchi: bool = False
    msmolbar: bool = False
    mskeepdir: bool = False
    msfragdist: WP = 0.0
    
    bhess: bool = True
    notemp: bool = True
    hotip: bool = False
    fermi: bool = False

    nots: bool = False
    esim: bool = False
    plot: bool = False
    oplot: bool = False
    int_masses: bool = False
    noiso: bool = False # do not plot isotope pattern
    calcKER: bool = False
    
    # special modes
    cneintscale: bool = False
    dxtb: bool = False

    sortoutcascade: bool = False
    ignoreip: bool = False
    reoptts: bool = True
    tsoptact: bool = False
    checkmult: bool = False
    pathmult: bool = False
    eyring: bool = True
    eyzpve: bool = False
    mthr: WP = 0.0
    intthr: WP = 0.0
    sthr: int = 150
    ithr: int = 100
    
    nfrag: int = 3 # number of subsequent fragmentations, at max 3
    
    # special stuff
    noirc: bool = False
    printlevel: int = 1
    logo: bool = False
    removedirs: bool = True
    
    # for timings (will be handled by Timer class instances, not directly here)
    tcrest: WP = 0.0 
    tip: WP = 0.0
    tpath: WP = 0.0
    tts: WP = 0.0
    thess: WP = 0.0
    tbhess: WP = 0.0
    
    # for testing
    tfscale: WP = 0.0
    nthermosteps: int = 200
    noplotisos: bool = False
    ircrun: bool = False
    picktsrun: bool = False

    # CID DATA
    prot: bool = False
    deprot: bool = False
    cid_mode: int = 1
    cid_elab: WP = 40.0
    cid_esi: WP = 0.0
    cid_esiatom: WP = 0.0
    cid_esiw: WP = 0.2
    cid_collw: WP = 0.5
    cid_maxcoll: int = 10
    cid_lchamb: WP = 0.25
    cid_pgas: WP = 0.132 # formerly 0.000132
    cid_TGas: WP = 300.0
    cid_mgas: WP = 39.94800 # argon
    cid_rgas: WP = 3.55266638 # bohr vdw-radius of Ar
    cid_scool: WP = 1.0
    solv: bool = False

# Utility functions
def eval_time(total_seconds: WP) -> str:
    """Transform a time in seconds to a string in the format hh:mm:ss"""
    hours = int(total_seconds // 3600)
    remaining_seconds = total_seconds % 3600
    minutes = int(remaining_seconds // 60)
    seconds = int(remaining_seconds % 60)
    return f"{hours}h :{minutes:02}m :{seconds:02}s"

def eval_sub_timer(timer_instance: Timer):
    if not timer_instance or timer_instance.times < 1:
        return
    for i in range(timer_instance.times):
        accumulated_time = timer_instance.t[i][2]
        if accumulated_time > 0:
            formatted_time = eval_time(accumulated_time)
            if timer_instance.names[i]:
                print(f"{timer_instance.names[i]:<24} wall time : {formatted_time:<20}")
        else:
            continue # No time recorded or name is empty

def eval_timer(timer_instance: Timer):
    print("\nWall Time Summary")
    small_head("Wall Time Summary") # Corrected to use the new function name
    eval_sub_timer(timer_instance)
    
    total_accumulated_time = sum(timer_instance.t[i][2] for i in range(timer_instance.times) if timer_instance.t)
    formatted_total_time = eval_time(total_accumulated_time)
    
    print("--------------------")
    print(f"Overall wall time  : {formatted_total_time}")
    
    timer_instance.clear()

def prop_quit(timer_instance: Timer): # Corrected to use the new function name
    eval_timer(timer_instance)
    print("\nQCxMS2 terminated normally.")
    sys.exit()

def small_head(text: str): # Corrected to use the new function name
    strlen = len(text)
    line = '-' * strlen
    print(f"\n {line}")
    print(f" {text}")
    print(f" {line}")

# Example of how qm_set%tcrest would be handled:
# Create a global timer instance or pass it around.
# main_timer = Timer()
# main_timer.init(10) # example: 10 different timings to track
# ...
# main_timer.start(0, "crest_calculation")
# ... call CREST ...
# main_timer.stop(0)
# run_type_data_instance.tcrest = main_timer.t[0][2] # Store accumulated time if needed
# ...
# prop_quit(main_timer)

```
