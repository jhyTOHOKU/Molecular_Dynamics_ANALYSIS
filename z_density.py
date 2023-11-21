from collections import defaultdict
import numpy as np
import pathlib

def process_trajectory_file_for_z_density(filepath, water_type, radius, z_min, z_max):
    positions_z = []
    with open(filepath, 'r') as file:
        for line in file:
            tokens = line.split()
            if len(tokens) == 9 and int(tokens[2]) == water_type:
                x, y, z = map(float, tokens[3:6])
                if z_min <= z < z_max and x < 0 and x**2 + y**2 <= radius**2:
                    positions_z.append(z)
    return positions_z

def calculate_z_density(all_positions_z, z_min, z_max, bin_width, molar_mass_water, avogadro_number, area):
    num_bins = int((z_max - z_min) / bin_width)
    density_dict = defaultdict(int)

    for z in all_positions_z:
        bin_index = int((z - z_min) / bin_width)
        if 0 <= bin_index < num_bins:
            density_dict[bin_index] += 1

    z_density = {}
    for bin_index in range(num_bins):
        bin_center = z_min + bin_width * (bin_index + 0.5)
        count = density_dict[bin_index]
        volume = bin_width * area
        number_density = count / volume
        mass_density = number_density * molar_mass_water / avogadro_number * 1e30
        z_density[bin_center] = mass_density

    return z_density

def calculate_and_write_z_density(filepath, output_dir, water_type, radius, z_min, z_max, bin_width, molar_mass_water, avogadro_number, bulk_water_density):
    positions_z = process_trajectory_file_for_z_density(filepath, water_type, radius, z_min, z_max)
    area = np.pi * radius**2 / 2
    z_density = calculate_z_density(positions_z, z_min, z_max, bin_width, molar_mass_water, avogadro_number, area)

    # Create output directory if it doesn't exist
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Derive output file name from trajectory file
    output_filename = pathlib.Path(output_dir) / (pathlib.Path(filepath).stem + '_z_density_normalized.csv')

    with open(output_filename, 'w') as outfile:
        outfile.write("Z_(Angstrom),Normalized_mass_density_(g/m^3)\n")
        # Sort the bin centers by their absolute values in ascending order
        for bin_center in sorted(z_density, key=lambda x: abs(x)):
            normalized_density = z_density[bin_center] / bulk_water_density
            # Use the absolute value of the bin center for output
            outfile.write(f"{abs(bin_center):.5f},{normalized_density:.5f}\n")

def calculate_from_files(trajectory_dir, output_dir, start_time, end_time, interval, water_type, radius, z_min, z_max, bin_width, molar_mass_water, avogadro_number, bulk_water_density):
    all_files = sorted(pathlib.Path(trajectory_dir).glob('*.lammpstrj'))
    files_to_analyze = [f for f in all_files if start_time <= int(f.stem.split('_')[1]) <= end_time and (int(f.stem.split('_')[1]) - start_time) % interval == 0]
    
    for file_number, filepath in enumerate(files_to_analyze, start=1):
        calculate_and_write_z_density(str(filepath), output_dir, water_type, radius, z_min, z_max, bin_width, molar_mass_water, avogadro_number, bulk_water_density)
        progress = (file_number / len(files_to_analyze)) * 100
        print(f"Processed file {file_number}/{len(files_to_analyze)}: {filepath.name} - Progress: {progress:.2f}%")


calculate_from_files(
    trajectory_dir='./trajectory',
    output_dir='./output_analysis/z_density',
    start_time=39999000,
    end_time=40000000,
    interval=1000,
    water_type=3,
    radius=13.560057102018412,
    z_min=-75,
    z_max=0,
    bin_width=3,
    molar_mass_water=18.01528,
    avogadro_number=6.02214076e23,
    bulk_water_density=10**6
)
