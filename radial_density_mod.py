from collections import defaultdict
import numpy as np
import pathlib

def process_trajectory_file(filepath, water_type, cutoff_radius):
    with open(filepath, 'r') as file:
        positions = []
        for line in file:
            tokens = line.split()
            if len(tokens) == 9 and int(tokens[2]) == water_type:
                x, y, z = map(float, tokens[3:6])
                if -80 < z < 0:
                    r = (x**2 + y**2)**0.5
                    if r <= cutoff_radius:
                        positions.append(r)
    return positions

def calculate_radial_density(all_positions, cutoff_radius, bin_width, num_files_to_analyze, molar_mass_water, avogadro_number, bulk_water_density):
    num_bins = int(cutoff_radius / bin_width) + 1
    density_dict = defaultdict(int)

    for r in all_positions:
        bin_index = int(r / bin_width)
        if bin_index < num_bins:
            density_dict[bin_index * bin_width] += 1

    average_density = {}
    for bin_center in np.arange(0, cutoff_radius, bin_width):
        count = density_dict[bin_center]
        volume = np.pi * ((bin_center + bin_width)**2 - bin_center**2) * 80
        number_density = count / volume / num_files_to_analyze
        mass_density = number_density * molar_mass_water / avogadro_number * 1e30
        normalized_density = mass_density / bulk_water_density
        average_density[bin_center + bin_width / 2] = normalized_density

    return average_density

def calculate_from_files(trajectory_dir, output_dir, output_filename, start_time, end_time, interval, water_type, cutoff_radius, bin_width, molar_mass_water, avogadro_number, bulk_water_density):
    all_files = sorted(pathlib.Path(trajectory_dir).glob('*.lammpstrj'))
    files_to_analyze = [f for f in all_files if start_time <= int(f.stem.split('_')[1]) <= end_time and (int(f.stem.split('_')[1]) - start_time) % interval == 0]
    num_files_to_analyze = len(files_to_analyze)

    all_positions = []
    for file_number, filepath in enumerate(files_to_analyze, start=1):
        all_positions.extend(process_trajectory_file(filepath, water_type, cutoff_radius))
        progress = (file_number / num_files_to_analyze) * 100
        print(f"Processed file {file_number}/{num_files_to_analyze}: {filepath.name} - Progress: {progress:.2f}%", flush=True)

    average_density = calculate_radial_density(all_positions, cutoff_radius, bin_width, num_files_to_analyze, molar_mass_water, avogadro_number, bulk_water_density)
    
    with open(pathlib.Path(output_dir) / output_filename, 'w') as outfile:
        outfile.write("Radial_distance_(Angstrom),Normalized_mass_density\n")
        for bin_center, normalized_density in average_density.items():
            # outfile.write(f"{bin_center:.2f},{normalized_density:.5f}\n")
            outfile.write(f"{bin_center:.2f},{normalized_density:.5f}\n")
