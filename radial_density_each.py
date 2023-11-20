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

def calculate_radial_density(all_positions, cutoff_radius, bin_width, molar_mass_water, avogadro_number, bulk_water_density):
    num_bins = int(cutoff_radius / bin_width) + 1
    density_dict = defaultdict(int)

    for r in all_positions:
        bin_index = int(r / bin_width)
        if bin_index < num_bins:
            density_dict[bin_index * bin_width] += 1

    average_density = {}
    for bin_center in np.arange(0, cutoff_radius, bin_width):
        count = density_dict[bin_center]
        volume = np.pi * ((bin_center + bin_width)**2 - bin_center**2) * (1.42000+73.97560)  
        number_density = count / volume
        mass_density = number_density * molar_mass_water / avogadro_number * 1e30
        normalized_density = mass_density / bulk_water_density
        average_density[bin_center + bin_width / 2] = normalized_density

    return average_density


def calculate_and_write_radial_density(filepath, output_dir, water_type, cutoff_radius, bin_width, molar_mass_water, avogadro_number, bulk_water_density):
    positions = process_trajectory_file(filepath, water_type, cutoff_radius)
    density = calculate_radial_density(positions, cutoff_radius, bin_width, molar_mass_water, avogadro_number, bulk_water_density)

    # Create output directory if it doesn't exist
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Derive output file name from trajectory file
    output_filename = pathlib.Path(output_dir) / (pathlib.Path(filepath).stem + '_radial_density.csv')
    
    with open(output_filename, 'w') as outfile:
        outfile.write("Radial_distance_(Angstrom),Normalized_mass_density\n")
        for bin_center, normalized_density in density.items():
            outfile.write(f"{bin_center:.5f},{normalized_density:.5f}\n")


def calculate_from_files(trajectory_dir, output_dir, start_time, end_time, interval, water_type, cutoff_radius, bin_width, molar_mass_water, avogadro_number, bulk_water_density):
    all_files = sorted(pathlib.Path(trajectory_dir).glob('*.lammpstrj'))
    files_to_analyze = [f for f in all_files if start_time <= int(f.stem.split('_')[1]) <= end_time and (int(f.stem.split('_')[1]) - start_time) % interval == 0]

    for file_number, filepath in enumerate(files_to_analyze, start=1):
        calculate_and_write_radial_density(str(filepath), output_dir, water_type, cutoff_radius, bin_width, molar_mass_water, avogadro_number, bulk_water_density)
        progress = (file_number / len(files_to_analyze)) * 100
        print(f"Processed file {file_number}/{len(files_to_analyze)}: {filepath.name} - Progress: {progress:.2f}%", flush=True)

#below : alterable 

calculate_from_files(
    trajectory_dir = './trajectory',
    output_dir = './output_analysis/water_radial_density',
    start_time = 39000000,
    end_time = 40000000,
    interval = 1000,
    water_type = 3,  # Assuming oxygen(water) type is 3
    cutoff_radius = 30,
    bin_width = 0.05,
    molar_mass_water = 18.01528,
    avogadro_number = 6.02214076e23,
    bulk_water_density = 1e6
)
