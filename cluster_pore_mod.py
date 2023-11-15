import MDAnalysis as mda
from MDAnalysis.analysis import distances
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import ConvexHull
import numpy as np
import os
import csv

def calculate_volume_and_area(cluster):
    if len(cluster) < 4:
        return np.nan, np.nan  # Not enough points to form a hull
    try:
        hull = ConvexHull(cluster)
        return hull.volume, hull.area
    except Exception:
        return np.nan, np.nan

def main():
    u = mda.Universe('system.data')  # Ensure this file is accessible

    start_time, end_time = 18e6, 20e6
    time_step = 1000  # Adjust as needed

    water = u.select_atoms("type 3 4")       # Update to match your topology
    hydronium = u.select_atoms("type 1 2")   # Update to match your topology
    all_molecules = water + hydronium

    output_dir = "./output_analysis/cluster_analysis_pore"
    os.makedirs(output_dir, exist_ok=True)

    with open('progress.log', 'w') as log_file:
        for time in range(int(start_time), int(end_time) + time_step, time_step):
            traj_file = f"trajectory/traj_{time:010d}.lammpstrj"
            u.load_new(traj_file, format='LAMMPSDUMP')
            log_file.write(f"Processing {traj_file}\n")
            log_file.flush()

            cluster_data = []
            for ts in u.trajectory:
                if start_time <= ts.time <= end_time:
                    distance_matrix = distances.distance_array(all_molecules.positions, all_molecules.positions, box=u.dimensions)
                    condensed_matrix = squareform(distance_matrix, force='tovector')
                    Z = linkage(condensed_matrix, method='single')
                    labels = fcluster(Z, t=3.5, criterion='distance')

                    for cluster_id in np.unique(labels):
                        cluster_atoms = all_molecules[labels == cluster_id]
                        max_z = np.max(cluster_atoms.positions[:, 2]) - 80  # Adjusted to simulation system range
                        min_z = np.min(cluster_atoms.positions[:, 2]) - 80  # Adjusted to simulation system range

                        if min_z > 0:
                            continue

                        volume, area = calculate_volume_and_area(cluster_atoms.positions)
                        cluster_data.append([len(cluster_atoms), volume, area, max_z, min_z])

            # Sort clusters by size (largest first)
            cluster_data.sort(key=lambda x: x[0], reverse=True)

            # Write to CSV file
            output_file = os.path.join(output_dir, f"cluster_analysis_{time:010d}.csv")
            with open(output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["size(atoms)","Volume(10nm^3)", "surface_area(10nm^2)", "max_z(10nm)", "min_z(10nm)"])
                writer.writerows(cluster_data)

            log_file.write(f"Completed processing for {traj_file}\n")
            log_file.flush()

if __name__ == "__main__":
    main()

#19998000