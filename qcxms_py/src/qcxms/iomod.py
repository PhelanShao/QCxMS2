import os
import shutil
import subprocess
from pathlib import Path
from typing import TypeVar, Union, List, Optional, Tuple

# Attempt to import WP, otherwise default to float
try:
    from .data import WP
except ImportError:
    WP = float

# --- File System Operations ---

def remove_file(filepath: Union[str, Path]):
    """Deletes a file. If filepath contains '*', the operation is skipped."""
    if isinstance(filepath, Path):
        filepath_str = str(filepath)
    else:
        filepath_str = filepath

    if "*" in filepath_str:
        # print(f"Skipping removal of wildcard path: {filepath_str}")
        return
    try:
        Path(filepath_str).unlink(missing_ok=True) # missing_ok=True mimics Fortran's lack of error if not found after open/close
    except IsADirectoryError:
        print(f"Warning: Attempted to remove directory {filepath_str} with remove_file. Use rmrf instead.")
    except OSError as e:
        print(f"Error removing file {filepath_str}: {e}")


def rmrf(path: Union[str, Path]):
    """Removes a file or directory recursively. Equivalent to 'rm -rf'."""
    path_obj = Path(path)
    try:
        if path_obj.is_file() or path_obj.is_symlink():
            path_obj.unlink(missing_ok=True)
        elif path_obj.is_dir():
            shutil.rmtree(path_obj, ignore_errors=True)
    except Exception as e: # Catch any potential error during deletion
        print(f"Error in rmrf for {path}: {e}")


def rmrf_wildcard(dataname_base: str):
    """Removes files or directories matching 'dataname_base*'."""
    # This is trickier due to shell globbing vs Python globbing.
    # 'rm -rf dataname_base*' in shell is powerful.
    # Path.glob is safer.
    # For simplicity, we'll execute the shell command if wildcards are essential.
    # Otherwise, one might list and then rmrf.
    # The Fortran seems to imply simple suffix wildcard.
    # This is a potential area for different behavior if not careful.
    # Let's assume it means to remove all files/dirs starting with dataname_base in the CWD.
    
    # Safer Pythonic way:
    # for p in Path(".").glob(f"{Path(dataname_base).name}*"): # Use name part for glob in current dir
    #     rmrf(p)
    # Fortran version likely just called `rm -rf dataname_base*` via system call.
    # Replicating with subprocess for now, assuming simple wildcard usage.
    try:
        # Warning: shell=True can be risky if dataname_base is untrusted.
        subprocess.run(f"rm -rf {dataname_base}*", shell=True, check=False)
    except Exception as e:
        print(f"Error in rmrf_wildcard for {dataname_base}*: {e}")


def touch(filepath: Union[str, Path]):
    """Creates an empty file (like touch) or updates its timestamp."""
    try:
        Path(filepath).touch()
    except OSError as e:
        print(f"Error touching file {filepath}: {e}")

def copy_file(source: Union[str, Path], destination: Union[str, Path]):
    """Copies a file from source to destination."""
    try:
        shutil.copy2(source, destination) # copy2 preserves more metadata
    except FileNotFoundError:
        # print(f"Warning: Source file {source} not found for copy.")
        pass # Fortran version also checks exist then copies, so no error if not found
    except shutil.SameFileError:
        # print(f"Warning: Source and destination are the same file: {source}")
        pass
    except OSError as e:
        print(f"Error copying file from {source} to {destination}: {e}")


def copy_to_subdir(source_file: Union[str, Path], dest_subdir: Union[str, Path]):
    """Copies a file into a subdirectory, keeping the original filename."""
    source_path = Path(source_file)
    dest_dir_path = Path(dest_subdir)
    
    if not source_path.exists():
        # print(f"Warning: Source file {source_file} not found for copy_to_subdir.")
        return

    if not dest_dir_path.is_dir():
        # print(f"Warning: Destination subdirectory {dest_subdir} does not exist or is not a directory.")
        # Fortran version changes dir, then opens file by original name.
        # This implies dest_subdir must exist.
        return

    try:
        shutil.copy2(source_path, dest_dir_path / source_path.name)
    except Exception as e:
        print(f"Error in copy_to_subdir from {source_file} to {dest_subdir}: {e}")


def move_file(source: Union[str, Path], destination: Union[str, Path]):
    """Moves/renames a file. Destination is overwritten if it exists."""
    source_path = Path(source)
    dest_path = Path(destination)
    try:
        if not source_path.exists():
            # print(f"File {source_path} does not exist for move!") # Matched Fortran output
            return
        if dest_path.exists():
            dest_path.unlink(missing_ok=True) # Try to remove destination first
        source_path.rename(dest_path)
    except OSError as e:
        print(f"Error moving file from {source} to {destination}: {e}")


def get_absolute_path(relative_path: Union[str, Path]) -> str:
    """Converts a relative path to an absolute path."""
    # Fortran version temporarily changed directory. Python's Path.resolve() is better.
    try:
        return str(Path(relative_path).resolve(strict=False)) # strict=False allows non-existent paths
    except Exception as e: # Catch potential errors during resolution
        print(f"Error resolving path {relative_path}: {e}")
        return str(Path(relative_path)) # Return original path string on error

def inquire_in_subdir(filename: str, subdir: Union[str, Path]) -> bool:
    """Checks if a file exists in a given subdirectory."""
    return (Path(subdir) / filename).exists()

def make_directory(dir_path: Union[str, Path], mode: int = 0o770):
    """Creates a directory. Mode is in octal (e.g., 0o770)."""
    try:
        Path(dir_path).mkdir(mode=mode, parents=True, exist_ok=True)
        return 0 # Success
    except OSError as e:
        print(f"Error creating directory {dir_path}: {e}")
        return -1 # Failure, similar to C mkdir return

def create_symlink(source: Union[str, Path], link_name: Union[str, Path]):
    """Creates a symbolic link pointing to source."""
    try:
        os.symlink(source, link_name)
        return 0 # Success
    except OSError as e:
        # print(f"Error creating symlink from {source} to {link_name}: {e}")
        return -1 # Failure

# --- File I/O & Parsing ---

def cat_file(filepath: Union[str, Path]):
    """Prints the content of a file to standard output."""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                print(line.rstrip('\n')) # rstrip to match Fortran's trim behavior on print
    except FileNotFoundError:
        print(f"Error: File {filepath} not found for cat_file.")
    except Exception as e:
        print(f"Error reading file {filepath} in cat_file: {e}")

def catdel_file(filepath: Union[str, Path]):
    """Prints file content then deletes the file."""
    cat_file(filepath)
    remove_file(filepath)


def append_to_file(source_file: Union[str, Path], target_file: Union[str, Path]):
    """Appends content of source_file to target_file."""
    try:
        with open(source_file, 'r') as f_source, open(target_file, 'a') as f_target:
            for line in f_source:
                f_target.write(line)
    except FileNotFoundError:
        print(f"Error: One or both files not found for append_to_file ({source_file}, {target_file}).")
    except Exception as e:
        print(f"Error during append_to_file: {e}")


def wrshort_int(filepath: Union[str, Path], value: int):
    """Writes a single integer to a file, overwriting it."""
    try:
        with open(filepath, 'w') as f:
            f.write(f"{value}\n")
    except IOError as e:
        print(f"Error writing integer to {filepath}: {e}")

def wrshort_real(filepath: Union[str, Path], value: WP):
    """Writes a single float (WP) to a file, overwriting it."""
    try:
        with open(filepath, 'w') as f:
            f.write(f"{value}\n") # Standard float formatting
    except IOError as e:
        print(f"Error writing real to {filepath}: {e}")

def wrshort_string(filepath: Union[str, Path], text: str):
    """Writes a single string to a file, overwriting it."""
    try:
        with open(filepath, 'w') as f:
            f.write(f"{text}\n")
    except IOError as e:
        print(f"Error writing string to {filepath}: {e}")


def rdshort_int(filepath: Union[str, Path], default: int = 0) -> int:
    """Reads a single integer from the first line of a file."""
    try:
        with open(filepath, 'r') as f:
            line = f.readline()
            return int(line.strip())
    except (FileNotFoundError, ValueError, IOError):
        return default

def rdshort_real(filepath: Union[str, Path], default: WP = 0.0) -> WP:
    """Reads a single float (WP) from the first line of a file."""
    try:
        with open(filepath, 'r') as f:
            line = f.readline()
            return float(line.strip())
    except (FileNotFoundError, ValueError, IOError):
        return default

def rdshort_string(filepath: Union[str, Path], default: str = "") -> str:
    """Reads a single string (first line) from a file."""
    try:
        with open(filepath, 'r') as f:
            return f.readline().rstrip('\n') # rstrip to remove newline, keep other whitespace
    except (FileNotFoundError, IOError):
        return default

def rdshort2_real(filepath: Union[str, Path], default: WP = 0.0) -> WP:
    """Reads a single float (WP) from the *second* line of a file."""
    try:
        with open(filepath, 'r') as f:
            f.readline() # Skip first line
            line = f.readline()
            return float(line.strip())
    except (FileNotFoundError, ValueError, IOError, AttributeError): # AttributeError if readline returns None
        return default

def minigrep(filepath: Union[str, Path], search_string: str) -> bool:
    """Checks if search_string is present in the file."""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if search_string in line:
                    return True
        return False
    except (FileNotFoundError, IOError):
        return False

def grepval(filepath: Union[str, Path], search_string: str, default_val: WP = 0.0) -> Tuple[bool, WP]:
    """
    Finds search_string in file. If found, returns True and the float
    value immediately following it on the same line. Otherwise, False and default_val.
    """
    try:
        with open(filepath, 'r') as f:
            for line in f:
                idx = line.find(search_string)
                if idx != -1:
                    # Extract substring after search_string
                    val_str_part = line[idx + len(search_string):].strip()
                    # Take the first whitespace-separated token from this part
                    val_token = val_str_part.split(maxsplit=1)[0] if val_str_part else ""
                    try:
                        return True, float(val_token)
                    except ValueError:
                        # Found string, but no valid float followed. Fortran might loop or error.
                        # For now, first match attempt. If it's not a float, we fail for this line.
                        continue # Try next line if first part wasn't a float
            return False, default_val # String not found or no valid float after it on any line
    except (FileNotFoundError, IOError):
        return False, default_val


def glinex(line: str, x: int) -> str:
    """
    Returns the x-th word (1-indexed) from a line.
    If x < 0, returns the last word.
    """
    words = line.split()
    if not words:
        return ""
    
    if x < 0: # Get last element
        return words[-1]
    if 1 <= x <= len(words):
        return words[x - 1] # Adjust to 0-indexed
    return "" # x out of bounds

def clinex(line: str, x: int) -> str:
    """Cuts the leading x words (1-indexed) from a string and returns the rest."""
    words = line.split()
    if x <= 0:
        return line # No words to cut or invalid x
    if x >= len(words):
        return "" # Cut all words or more
    
    return " ".join(words[x:])


def clear_setblock_in_file(filepath: Union[str, Path]):
    """
    Reads a file, copies content to a temp file until '$set' or '$end' is encountered,
    writes '$end', then replaces the original file.
    """
    temp_filepath = Path(str(filepath) + ".setdgtmp")
    original_path = Path(filepath)

    if not original_path.exists():
        # print(f"File {filepath} not found for clear_setblock.")
        return

    try:
        with open(original_path, 'r') as f_in, open(temp_filepath, 'w') as f_out:
            for line_content in f_in:
                stripped_line = line_content.strip().lower() # For case-insensitive match
                if "$set" in stripped_line or "$end" in stripped_line:
                    f_out.write("$end\n")
                    break
                else:
                    f_out.write(line_content)
        
        # Replace original file with temp file
        shutil.move(str(temp_filepath), str(original_path))

    except IOError as e:
        print(f"Error during clear_setblock for {filepath}: {e}")
        if temp_filepath.exists():
            temp_filepath.unlink(missing_ok=True) # Clean up temp file on error
    except Exception as e: # Catch any other unexpected error
        print(f"Unexpected error during clear_setblock for {filepath}: {e}")
        if temp_filepath.exists():
            temp_filepath.unlink(missing_ok=True)


def file_checker(filepath_in: Union[str, Path]) -> Tuple[bool, str]:
    """Checks if filepath_in exists. Returns (exists_bool, path_str_or_empty)."""
    path_obj = Path(filepath_in)
    if path_obj.exists():
        return True, str(path_obj)
    else:
        return False, ""

# --- Environment Variables ---

def set_env_var(var_name: str, value: Union[str, int, float]) -> int:
    """Sets an environment variable. Returns 0 on success, -1 on error."""
    try:
        os.environ[var_name] = str(value)
        return 0
    except Exception as e:
        # print(f"Error setting environment variable {var_name}: {e}")
        return -1 # C setenv returns -1 on error

# --- String Manipulation ---

def to_upper(input_string: str) -> str:
    """Converts a string to uppercase."""
    return input_string.upper()

def to_lower(input_string: str) -> str:
    """Converts a string to lowercase."""
    return input_string.lower()

# --- Miscellaneous ---

def password_lock(password_reference: str) -> bool:
    """A simple password lock mechanism. Returns True if password matches, False otherwise."""
    # In a real application, use `getpass` module for passwords.
    # This is a direct translation of the Fortran behavior.
    try:
        entered_password = input(f"Locked feature. Enter password: ")
        if entered_password == password_reference:
            print("Valid. Continue.\n")
            return True
        else:
            # The Fortran code calls `error stop`. Here we'll return False and let caller decide.
            print("Invalid.\n") 
            return False
    except Exception: # Handle potential errors during input
        print("Error during password entry.\n")
        return False

# --- Centralized Command Execution ---
def execute_command(command_list: List[str], cwd: Optional[Union[str, Path]] = None, 
                    shell: bool = False, capture_output: bool = True, text: bool = True) -> subprocess.CompletedProcess:
    """
    A more robust version of execute_command_line from utility.
    command_list: A list of command arguments (e.g., ['ls', '-l'])
    """
    try:
        result = subprocess.run(
            command_list if not shell else " ".join(command_list), 
            cwd=cwd, 
            shell=shell, 
            capture_output=capture_output, 
            text=text,
            check=False # Don't raise exception for non-zero exit codes by default
        )
        return result
    except FileNotFoundError: # Command not found
        print(f"Error: Command '{command_list[0]}' not found.")
        return subprocess.CompletedProcess(args=command_list, returncode=127, stdout="", stderr="Command not found")
    except Exception as e:
        print(f"Error executing command '{' '.join(command_list)}': {e}")
        return subprocess.CompletedProcess(args=command_list, returncode=-1, stdout="", stderr=str(e))


if __name__ == '__main__':
    # Example usage:
    test_dir = Path("iomod_test_dir")
    test_dir.mkdir(exist_ok=True)

    # Test wrshort and rdshort
    wrshort_string(test_dir / "test.txt", "Hello World")
    print(f"rdshort_string: {rdshort_string(test_dir / 'test.txt')}")
    wrshort_int(test_dir / "test.int", 123)
    print(f"rdshort_int: {rdshort_int(test_dir / 'test.int')}")
    wrshort_real(test_dir / "test.real", 45.67)
    print(f"rdshort_real: {rdshort_real(test_dir / 'test.real')}")
    print(f"rdshort2_real (from test.real, should be 0): {rdshort2_real(test_dir / 'test.real')}")


    # Test grepval
    with open(test_dir / "grep_test.txt", "w") as f:
        f.write("Value is 123.45\n")
        f.write("Another value is: -0.56 \n")
        f.write("No value here\n")
        f.write("energy = 3.14159\n")
    
    found, val = grepval(test_dir / "grep_test.txt", "Value is ")
    print(f"grepval 'Value is ': Found={found}, Val={val}") # Exp: True, 123.45
    found, val = grepval(test_dir / "grep_test.txt", "value is: ") # case diff
    print(f"grepval 'value is: ': Found={found}, Val={val}") # Exp: True, -0.56
    found, val = grepval(test_dir / "grep_test.txt", "energy = ")
    print(f"grepval 'energy = ': Found={found}, Val={val}") # Exp: True, 3.14159
    found, val = grepval(test_dir / "grep_test.txt", "nonexistent ")
    print(f"grepval 'nonexistent ': Found={found}, Val={val}") # Exp: False, 0.0

    # Test glinex
    line_for_glinex = "  word1  word2   word3  word4 "
    print(f"glinex '{line_for_glinex}', 1: '{glinex(line_for_glinex, 1)}'") # word1
    print(f"glinex '{line_for_glinex}', 3: '{glinex(line_for_glinex, 3)}'") # word3
    print(f"glinex '{line_for_glinex}', -1: '{glinex(line_for_glinex, -1)}'")# word4
    print(f"glinex '{line_for_glinex}', 5: '{glinex(line_for_glinex, 5)}'") # ""

    # Test clinex
    print(f"clinex '{line_for_glinex}', 1: '{clinex(line_for_glinex, 1)}'") # "word2   word3  word4 "
    print(f"clinex '{line_for_glinex}', 2: '{clinex(line_for_glinex, 2)}'") # "word3  word4 "
    
    # Test file ops
    touch(test_dir / "touched.txt")
    print(f"touched.txt exists: {(test_dir / 'touched.txt').exists()}")
    copy_file(test_dir / "test.txt", test_dir / "test_copy.txt")
    print(f"test_copy.txt exists: {(test_dir / 'test_copy.txt').exists()}")
    move_file(test_dir / "test_copy.txt", test_dir / "test_moved.txt")
    print(f"test_moved.txt exists: {(test_dir / 'test_moved.txt').exists()}")
    print(f"test_copy.txt exists after move: {(test_dir / 'test_copy.txt').exists()}")
    
    # Test execute_command
    res = execute_command(["echo", "Hello from execute_command"])
    print(f"execute_command echo: stdout='{res.stdout.strip()}', stderr='{res.stderr.strip()}', retcode={res.returncode}")

    # Cleanup
    shutil.rmtree(test_dir)
    print(f"Cleaned up {test_dir}")
