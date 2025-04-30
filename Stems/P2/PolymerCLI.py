from polymerClasses import Macromolecule
from statistics import mean, stdev


def get_input(prompt, default, type_cast=str):
    """Get user input with default fallback and type casting."""
    response = input(f"{prompt} ({default})?: ").strip().lower()
    try:
        return type_cast(response) if response else default
    except ValueError:
        print(f"Invalid input. Using default: {default}")
        return default


def run_simulation():
    """
    Command-line interface for simulating freely jointed polymer chains of poly(ethylene).
    Collects user input for degree of polymerization (N) and number of molecules,
    runs simulation, and prints out statistics like:
    - Center of mass
    - Radius of gyration
    - End-to-end distance
    - Polydispersity index (PDI)
    """

    target_N = get_input("Degree of polymerization", 1000, int)
    num_mols = get_input("How many molecules", 50, int)

    polymers = []
    for _ in range(num_mols):
        pm = Macromolecule(target_N)
        pm.build_chain()
        polymers.append(pm)

    # Convert outputs
    com_x = [p.com.x * 1e9 for p in polymers]  # to nanometers
    com_y = [p.com.y * 1e9 for p in polymers]
    com_z = [p.com.z * 1e9 for p in polymers]

    e2e_um = [p.e2e * 1e6 for p in polymers]  # to micrometers
    rog_um = [p.rog * 1e6 for p in polymers]
    weights = [p.weight for p in polymers]

    # Output
    print(f"\nMetrics for {num_mols} molecules of degree of polymerization = {target_N}")
    print(f"Avg. Center of Mass (nm) = {mean(com_x):.3f}, {mean(com_y):.3f}, {mean(com_z):.3f}")
    print("End-to-end distance (μm):")
    print(f"\tAverage = {mean(e2e_um):.3f}")
    print(f"\tStd. Dev. = {stdev(e2e_um):.3f}" if len(e2e_um) > 1 else "\tStd. Dev. = 0.000")
    print("Radius of gyration (μm):")
    print(f"\tAverage = {mean(rog_um):.3f}")
    print(f"\tStd. Dev. = {stdev(rog_um):.3f}" if len(rog_um) > 1 else "\tStd. Dev. = 0.000")
    print(f"PDI = {sum(w ** 2 for w in weights) / sum(weights) / mean(weights):.2f}")



def main():
    print("\n🧪 Poly(ethylene) Freely Jointed Chain Simulator 🧪\n")
    again = True
    while again:
        run_simulation()
        again_input = input("\nWould you like to run another simulation? (Y/N): ").strip().lower()
        again = again_input in ["y", "yes", "true"]


if __name__ == "__main__":
    main()
