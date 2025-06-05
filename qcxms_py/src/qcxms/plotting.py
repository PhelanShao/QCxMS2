import os
import math
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.backends.backend_pdf import PdfPages
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available. Plotting functionality will be limited.")

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility
    from .constants import AUTOEV, EVTOKCAL
except ImportError:
    # Fallbacks for standalone/testing
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    AUTOEV = 27.211386245988
    EVTOKCAL = 23.060547830618307


def plot_mass_spectrum_py(
    env: RunTypeData,
    output_dir: Path,
    spectrum_file: str = "allspec.dat",
    title: str = "QCxMS2 Mass Spectrum"
) -> bool:
    """
    Plot mass spectrum from QCxMS2 simulation results.
    
    Args:
        env: Runtime environment data
        output_dir: Directory containing spectrum data
        spectrum_file: Name of spectrum data file
        title: Plot title
    
    Returns:
        True if successful, False otherwise
    """
    if not MATPLOTLIB_AVAILABLE:
        print("  Error: matplotlib not available for mass spectrum plotting")
        return False
    
    try:
        spectrum_path = output_dir / spectrum_file
        if not spectrum_path.exists():
            print(f"  Error: Spectrum file {spectrum_path} not found")
            return False
        
        # Read spectrum data
        masses = []
        intensities = []
        
        with open(spectrum_path, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mass = float(parts[0])
                            intensity = float(parts[1])
                            masses.append(mass)
                            intensities.append(intensity)
                        except ValueError:
                            continue
        
        if not masses:
            print(f"  Error: No valid data found in {spectrum_path}")
            return False
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot spectrum as stem plot
        ax.stem(masses, intensities, basefmt=' ', linefmt='b-', markerfmt='bo')
        
        # Formatting
        ax.set_xlabel('m/z', fontsize=12)
        ax.set_ylabel('Relative Intensity (%)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Set y-axis to start from 0
        ax.set_ylim(bottom=0)
        
        # Add some padding to x-axis
        if masses:
            x_range = max(masses) - min(masses)
            ax.set_xlim(min(masses) - 0.05 * x_range, max(masses) + 0.05 * x_range)
        
        # Add annotation for molecular ion peak
        if masses and intensities:
            max_mass_idx = masses.index(max(masses))
            molecular_ion_mass = masses[max_mass_idx]
            ax.annotate(f'M+ = {molecular_ion_mass:.1f}', 
                       xy=(molecular_ion_mass, intensities[max_mass_idx]),
                       xytext=(10, 10), textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        
        plt.tight_layout()
        
        # Save plot
        output_file = output_dir / "mass_spectrum.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        # Also save as PDF
        pdf_file = output_dir / "mass_spectrum.pdf"
        plt.savefig(pdf_file, bbox_inches='tight')
        
        plt.close()
        
        print(f"  Mass spectrum plotted and saved to {output_file}")
        return True
        
    except Exception as e:
        print(f"  Error plotting mass spectrum: {e}")
        return False


def plot_energy_distribution_py(
    env: RunTypeData,
    output_dir: Path,
    energy_file: str = "iee_distribution.dat",
    title: str = "Internal Energy Distribution"
) -> bool:
    """
    Plot internal energy distribution.
    
    Args:
        env: Runtime environment data
        output_dir: Directory containing energy data
        energy_file: Name of energy distribution file
        title: Plot title
    
    Returns:
        True if successful, False otherwise
    """
    if not MATPLOTLIB_AVAILABLE:
        print("  Error: matplotlib not available for energy distribution plotting")
        return False
    
    try:
        energy_path = output_dir / energy_file
        if not energy_path.exists():
            print(f"  Error: Energy file {energy_path} not found")
            return False
        
        # Read energy distribution data
        energies = []
        probabilities = []
        
        with open(energy_path, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            energy = float(parts[0])
                            prob = float(parts[1])
                            energies.append(energy)
                            probabilities.append(prob)
                        except ValueError:
                            continue
        
        if not energies:
            print(f"  Error: No valid data found in {energy_path}")
            return False
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot as line and fill
        ax.plot(energies, probabilities, 'b-', linewidth=2, label='Probability Density')
        ax.fill_between(energies, probabilities, alpha=0.3, color='blue')
        
        # Formatting
        ax.set_xlabel('Internal Energy (eV)', fontsize=12)
        ax.set_ylabel('Probability Density', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add statistics
        if energies and probabilities:
            # Calculate mean energy
            mean_energy = sum(e * p for e, p in zip(energies, probabilities)) / sum(probabilities)
            ax.axvline(mean_energy, color='red', linestyle='--', linewidth=2, 
                      label=f'Mean = {mean_energy:.2f} eV')
            ax.legend()
        
        plt.tight_layout()
        
        # Save plot
        output_file = output_dir / "energy_distribution.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        pdf_file = output_dir / "energy_distribution.pdf"
        plt.savefig(pdf_file, bbox_inches='tight')
        
        plt.close()
        
        print(f"  Energy distribution plotted and saved to {output_file}")
        return True
        
    except Exception as e:
        print(f"  Error plotting energy distribution: {e}")
        return False


def plot_fragmentation_tree_py(
    env: RunTypeData,
    output_dir: Path,
    fragments_file: str = "allfrags.dat",
    title: str = "Fragmentation Tree"
) -> bool:
    """
    Plot fragmentation tree showing hierarchical fragmentation pathways.
    
    Args:
        env: Runtime environment data
        output_dir: Directory containing fragmentation data
        fragments_file: Name of fragments data file
        title: Plot title
    
    Returns:
        True if successful, False otherwise
    """
    if not MATPLOTLIB_AVAILABLE:
        print("  Error: matplotlib not available for fragmentation tree plotting")
        return False
    
    try:
        fragments_path = output_dir / fragments_file
        if not fragments_path.exists():
            print(f"  Error: Fragments file {fragments_path} not found")
            return False
        
        # Read fragmentation data
        fragments = []
        with open(fragments_path, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                for line in lines[1:]:  # Skip first line (count)
                    if line.strip():
                        fragments.append(line.strip())
        
        if not fragments:
            print(f"  Error: No fragments found in {fragments_path}")
            return False
        
        # Parse fragment hierarchy
        fragment_tree = _parse_fragment_hierarchy_py(fragments)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # Plot tree structure
        _plot_tree_structure_py(ax, fragment_tree, title)
        
        plt.tight_layout()
        
        # Save plot
        output_file = output_dir / "fragmentation_tree.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        pdf_file = output_dir / "fragmentation_tree.pdf"
        plt.savefig(pdf_file, bbox_inches='tight')
        
        plt.close()
        
        print(f"  Fragmentation tree plotted and saved to {output_file}")
        return True
        
    except Exception as e:
        print(f"  Error plotting fragmentation tree: {e}")
        return False


def _parse_fragment_hierarchy_py(fragments: List[str]) -> Dict[str, Any]:
    """Parse fragment paths into hierarchical tree structure."""
    tree = {"name": "M+", "children": [], "level": 0}
    
    for fragment_path in fragments:
        # Parse fragment path (e.g., "p0f1", "p1f2p0", etc.)
        levels = []
        current = ""
        
        for char in fragment_path:
            if char in 'pf':
                if current:
                    levels.append(current)
                current = char
            else:
                current += char
        
        if current:
            levels.append(current)
        
        # Add to tree
        current_node = tree
        for i, level in enumerate(levels):
            # Find or create child node
            child_found = False
            for child in current_node["children"]:
                if child["name"] == level:
                    current_node = child
                    child_found = True
                    break
            
            if not child_found:
                new_child = {"name": level, "children": [], "level": i + 1}
                current_node["children"].append(new_child)
                current_node = new_child
    
    return tree


def _plot_tree_structure_py(ax, tree: Dict[str, Any], title: str):
    """Plot tree structure using matplotlib."""
    # Calculate positions for nodes
    positions = {}
    _calculate_positions_py(tree, positions, 0, 0, 1.0)
    
    # Plot nodes and connections
    for node_id, (x, y, node) in positions.items():
        # Plot node
        circle = patches.Circle((x, y), 0.05, facecolor='lightblue', 
                               edgecolor='black', linewidth=1)
        ax.add_patch(circle)
        
        # Add label
        ax.text(x, y - 0.1, node["name"], ha='center', va='top', fontsize=8)
        
        # Plot connections to children
        for child in node["children"]:
            child_id = id(child)
            if child_id in positions:
                child_x, child_y, _ = positions[child_id]
                ax.plot([x, child_x], [y, child_y], 'k-', linewidth=1, alpha=0.7)
    
    # Set axis properties
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-1.1, 0.1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)


def _calculate_positions_py(node: Dict[str, Any], positions: Dict[int, Tuple[float, float, Dict]], 
                           x: float, y: float, width: float):
    """Calculate positions for tree nodes."""
    node_id = id(node)
    positions[node_id] = (x, y, node)
    
    if node["children"]:
        child_width = width / len(node["children"])
        start_x = x - width / 2 + child_width / 2
        
        for i, child in enumerate(node["children"]):
            child_x = start_x + i * child_width
            child_y = y - 0.3  # Move down one level
            _calculate_positions_py(child, positions, child_x, child_y, child_width * 0.8)


def plot_reaction_profile_py(
    env: RunTypeData,
    output_dir: Path,
    reaction_dir: str,
    title: str = "Reaction Energy Profile"
) -> bool:
    """
    Plot reaction energy profile showing reactants, TS, and products.
    
    Args:
        env: Runtime environment data
        output_dir: Directory containing reaction data
        reaction_dir: Name of reaction directory
        title: Plot title
    
    Returns:
        True if successful, False otherwise
    """
    if not MATPLOTLIB_AVAILABLE:
        print("  Error: matplotlib not available for reaction profile plotting")
        return False
    
    try:
        reaction_path = output_dir / reaction_dir
        if not reaction_path.exists():
            print(f"  Error: Reaction directory {reaction_path} not found")
            return False
        
        # Read energy data
        energies = {}
        
        # Reactant energy
        reactant_qmdata = reaction_path / "qmdata"
        if reactant_qmdata.exists():
            # Parse reactant energy from qmdata
            with open(reactant_qmdata, 'r') as f:
                for line in f:
                    if "sp" in line:
                        parts = line.split()
                        if len(parts) >= 4:
                            try:
                                energies["reactant"] = float(parts[-1])
                                break
                            except ValueError:
                                continue
        
        # TS energy
        ts_dir = reaction_path / "ts"
        if ts_dir.exists():
            ts_qmdata = ts_dir / "qmdata"
            if ts_qmdata.exists():
                with open(ts_qmdata, 'r') as f:
                    for line in f:
                        if "sp" in line:
                            parts = line.split()
                            if len(parts) >= 4:
                                try:
                                    energies["ts"] = float(parts[-1])
                                    break
                                except ValueError:
                                    continue
        
        # Product energy (simplified - would need to read from product directories)
        # For now, estimate from reaction energy
        de_file = reaction_path / f"de_{env.tslevel}"
        if de_file.exists():
            de_ev = iomod.rdshort_real(de_file, default=0.0)
            if "reactant" in energies:
                energies["product"] = energies["reactant"] + de_ev / AUTOEV
        
        if len(energies) < 2:
            print(f"  Error: Insufficient energy data found in {reaction_path}")
            return False
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Convert to relative energies (eV)
        ref_energy = energies.get("reactant", 0.0)
        rel_energies = {k: (v - ref_energy) * AUTOEV for k, v in energies.items()}
        
        # Plot energy profile
        x_positions = {"reactant": 0, "ts": 1, "product": 2}
        x_labels = ["Reactant", "TS", "Product"]
        
        x_vals = []
        y_vals = []
        labels = []
        
        for state in ["reactant", "ts", "product"]:
            if state in rel_energies:
                x_vals.append(x_positions[state])
                y_vals.append(rel_energies[state])
                labels.append(state.capitalize())
        
        # Plot points and lines
        ax.plot(x_vals, y_vals, 'bo-', linewidth=2, markersize=8)
        
        # Add horizontal lines for each state
        for i, (x, y) in enumerate(zip(x_vals, y_vals)):
            ax.hlines(y, x - 0.2, x + 0.2, colors='blue', linewidth=3)
            ax.text(x, y + 0.1, f'{y:.2f} eV', ha='center', va='bottom', fontsize=10)
        
        # Formatting
        ax.set_xlabel('Reaction Coordinate', fontsize=12)
        ax.set_ylabel('Relative Energy (eV)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Set x-axis labels
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels)
        
        # Add barrier annotation
        if "reactant" in rel_energies and "ts" in rel_energies:
            barrier = rel_energies["ts"] - rel_energies["reactant"]
            ax.annotate(f'Ea = {barrier:.2f} eV', 
                       xy=(0.5, (rel_energies["reactant"] + rel_energies["ts"]) / 2),
                       xytext=(10, 10), textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        
        plt.tight_layout()
        
        # Save plot
        output_file = output_dir / f"{reaction_dir}_profile.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        pdf_file = output_dir / f"{reaction_dir}_profile.pdf"
        plt.savefig(pdf_file, bbox_inches='tight')
        
        plt.close()
        
        print(f"  Reaction profile plotted and saved to {output_file}")
        return True
        
    except Exception as e:
        print(f"  Error plotting reaction profile: {e}")
        return False


def generate_all_plots_py(env: RunTypeData, output_dir: Path) -> bool:
    """
    Generate all available plots for QCxMS2 results.
    
    Args:
        env: Runtime environment data
        output_dir: Directory containing simulation results
    
    Returns:
        True if at least one plot was generated successfully
    """
    print(f"\n--- Generating QCxMS2 Plots ---")
    
    if not MATPLOTLIB_AVAILABLE:
        print("  Error: matplotlib not available. Cannot generate plots.")
        return False
    
    success_count = 0
    
    # Plot mass spectrum
    if plot_mass_spectrum_py(env, output_dir):
        success_count += 1
    
    # Plot energy distribution
    if plot_energy_distribution_py(env, output_dir):
        success_count += 1
    
    # Plot fragmentation tree
    if plot_fragmentation_tree_py(env, output_dir):
        success_count += 1
    
    # Plot reaction profiles for all reactions
    for reaction_dir in output_dir.glob("p*"):
        if reaction_dir.is_dir():
            if plot_reaction_profile_py(env, output_dir, reaction_dir.name):
                success_count += 1
    
    # Generate summary report
    if success_count > 0:
        _generate_plot_summary_py(env, output_dir, success_count)
    
    print(f"--- Plot Generation Complete: {success_count} plots created ---")
    return success_count > 0


def _generate_plot_summary_py(env: RunTypeData, output_dir: Path, plot_count: int):
    """Generate HTML summary of all plots."""
    try:
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>QCxMS2 Results Summary</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f0f0; padding: 10px; border-radius: 5px; }}
        .plot-section {{ margin: 20px 0; }}
        .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
        .plot-item {{ border: 1px solid #ddd; padding: 10px; border-radius: 5px; }}
        img {{ max-width: 100%; height: auto; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>QCxMS2 Simulation Results</h1>
        <p>Generated {plot_count} plots from simulation data</p>
        <p>Simulation parameters: Level = {env.tslevel}, Temperature = {env.temp} K</p>
    </div>
    
    <div class="plot-section">
        <h2>Mass Spectrum</h2>
        <div class="plot-item">
            <img src="mass_spectrum.png" alt="Mass Spectrum">
        </div>
    </div>
    
    <div class="plot-section">
        <h2>Energy Distribution</h2>
        <div class="plot-item">
            <img src="energy_distribution.png" alt="Energy Distribution">
        </div>
    </div>
    
    <div class="plot-section">
        <h2>Fragmentation Tree</h2>
        <div class="plot-item">
            <img src="fragmentation_tree.png" alt="Fragmentation Tree">
        </div>
    </div>
    
    <div class="plot-section">
        <h2>Reaction Profiles</h2>
        <div class="plot-grid">
"""
        
        # Add reaction profile images
        for img_file in output_dir.glob("p*_profile.png"):
            html_content += f"""
            <div class="plot-item">
                <h3>{img_file.stem}</h3>
                <img src="{img_file.name}" alt="{img_file.stem}">
            </div>
"""
        
        html_content += """
        </div>
    </div>
</body>
</html>
"""
        
        summary_file = output_dir / "results_summary.html"
        with open(summary_file, 'w') as f:
            f.write(html_content)
        
        print(f"  HTML summary generated: {summary_file}")
        
    except Exception as e:
        print(f"  Error generating plot summary: {e}")


if __name__ == '__main__':
    print("Testing plotting module...")
    
    # Create test environment
    env_test = RunTypeData()
    env_test.tslevel = "gfn2"
    env_test.temp = 298
    
    # Create test data directory
    test_dir = Path("test_plots")
    test_dir.mkdir(exist_ok=True)
    
    try:
        # Create dummy spectrum data
        with open(test_dir / "allspec.dat", 'w') as f:
            f.write("# Mass spectrum data\n")
            for i in range(10, 101, 10):
                intensity = 100 * math.exp(-(i-50)**2/500)
                f.write(f"{i} {intensity:.2f}\n")
        
        # Create dummy energy distribution
        with open(test_dir / "iee_distribution.dat", 'w') as f:
            f.write("# Energy distribution\n")
            for i in range(0, 51):
                energy = i * 0.1
                prob = math.exp(-energy/2) * energy
                f.write(f"{energy:.2f} {prob:.6f}\n")
        
        # Test plotting functions
        if MATPLOTLIB_AVAILABLE:
            print("Testing mass spectrum plotting...")
            plot_mass_spectrum_py(env_test, test_dir)
            
            print("Testing energy distribution plotting...")
            plot_energy_distribution_py(env_test, test_dir)
            
            print("Testing complete plot generation...")
            generate_all_plots_py(env_test, test_dir)
        else:
            print("Matplotlib not available - skipping plot tests")
    
    finally:
        # Cleanup
        import shutil
        if test_dir.exists():
            shutil.rmtree(test_dir)
    
    print("Plotting module testing completed.")