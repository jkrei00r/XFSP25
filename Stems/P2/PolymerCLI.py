from polymerClasses import Macromolecule
from statistics import mean, stdev


def main():
    target_N = int(input("degree of polymerization (1000)?: ") or 1000)
    num_mols = int(input("How many molecules (50)?: ") or 50)

    polymers = []
    for _ in range(num_mols):
        pm = Macromolecule(target_N)
        pm.build_chain()
        polymers.append(pm)

    # Convert to micrometers and nanometers
    com_x = [p.com.x * 1e9 for p in polymers]
    com_y = [p.com.y * 1e9 for p in polymers]
    com_z = [p.com.z * 1e9 for p in polymers]

    e2e_um = [p.e2e * 1e6 for p in polymers]
    rog_um = [p.rog * 1e6 for p in polymers]
    weights = [p.weight for p in polymers]

    print(f"\nMetrics for {num_mols} molecules of degree of polymerization = {target_N}")
    print(f"Avg. Center of Mass (nm) = {mean(com_x):.3f}, {mean(com_y):.3f}, {mean(com_z):.3f}")
    print("End-to-end distance (μm):")
    print(f"\tAverage = {mean(e2e_um):.3f}")
    print(f"\tStd. Dev. = {stdev(e2e_um):.3f}" if len(e2e_um) > 1 else "\tStd. Dev. = 0.000")
    print("Radius of gyration (μm):")
    print(f"\tAverage = {mean(rog_um):.3f}")
    print(f"\tStd. Dev. = {stdev(rog_um):.3f}" if len(rog_um) > 1 else "\tStd. Dev. = 0.000")
    print(f"PDI = {sum(w ** 2 for w in weights) / sum(weights) / mean(weights):.2f}")


if __name__ == "__main__":
    main()
