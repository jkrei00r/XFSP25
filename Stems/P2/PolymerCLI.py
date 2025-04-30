# PolymerCLI.py

"""
Command-Line Interface (CLI) for simulating freely jointed chains of poly(ethylene).

This script allows users to simulate multiple polymer chains and analyze their structural
properties, including:

- Center of Mass (COM)
- End-to-End Distance (Re)
- Radius of Gyration (Rg)
- Polydispersity Index (PDI)

Simulation assumes:
--------------------
- Freely jointed chain model in 3D space.
- Random orientation between monomers with a fixed bond length (1.54 nm).
- Gaussian distribution for chain length variability (10% standard deviation).

To run:
-------
$ python PolymerCLI.py

"""

from polymerClasses import Macromolecule
from statistics import mean, stdev


def get_input(prompt, default, type_cast=str):
    """
    Gets user input with default fallback and type casting.

    Parameters:
    -----------
    prompt : str
        The prompt message to display.
    default : any
        Default value if no input is provided.
    type_cast : type
        Function to cast the input into (e.g., int, float).

    Returns:
    --------
    type_cast(input) or default
    """
    response = input(f"{prompt} ({default})?: ").strip().lower()
    try:
        return type_cast(response) if response else default
    except ValueError:
        print(f"Invalid input. Using default: {default}")
        return default


def run_simulation():
    """
    Handles simulation logic:
    - Prompts user for polymer size and quantity
    - Runs simulation
    - Calculates and prints metrics for each set of chains
    """
    print("\n🔄 Starting Simulation...")

    # User inputs
    target_N = get_input("Enter target degree of polymerization (N)", 1000, int)
    num_mols = get_input("Enter number of molecules to simulate", 50, int)

    # Create and simulate polymers
    polymers = []
    for _ in range(num_mols):
        pm = Macromolecule(target_N)
        pm.build_chain()
        polymers.append(pm)

    # Extract metrics (convert units where needed)
    com_x = [p.com.x * 1e9 for p in polymers]  # nm
    com_y = [p.com.y * 1e9 for p in polymers]
    com_z = [p.com.z * 1e9 for p in polymers]

    e2e_um = [p.e2e * 1e6 for p in polymers]  # μm
    rog_um = [p.rog * 1e6 for p in polymers]
    weights = [p.weight for p in polymers]

    # Display results
    print(f"\n📊 Simulation Results for {num_mols} molecules with target N = {target_N}\n")
    print(f"🔸 Average Center of Mass (nm): X = {mean(com_x):.2f}, Y = {mean(com_y):.2f}, Z = {mean(com_z):.2f}")

    print("🔸 End-to-End Distance (μm):")
    print(f"\tMean = {mean(e2e_um):.4f}")
    print(f"\tStd Dev = {stdev(e2e_um):.4f}" if len(e2e_um) > 1 else "\tStd Dev = 0.0000")

    print("🔸 Radius of Gyration (μm):")
    print(f"\tMean = {mean(rog_um):.4f}")
    print(f"\tStd Dev = {stdev(rog_um):.4f}" if len(rog_um) > 1 else "\tStd Dev = 0.0000")

    # PDI calculation: Mw / Mn / Mavg
    pdi = sum(w ** 2 for w in weights) / sum(weights) / mean(weights)
    print(f"🔸 Polydispersity Index (PDI) = {pdi:.3f}")


def main():
    """
    Entry point for the CLI application. Repeats simulation until user opts out.
    """
    print("\n🧪 Poly(ethylene) Freely Jointed Chain Simulator 🧪")

    again = True
    while again:
        run_simulation()
        again_input = input("\nRun another simulation? (Y/N): ").strip().lower()
        again = again_input in ["y", "yes", "true"]


# Only run if executed as a script
if __name__ == "__main__":
    main()
