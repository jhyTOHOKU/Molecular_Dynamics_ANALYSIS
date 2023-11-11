import configparser
import sys
sys.path.append("intput path")
import radial_density_mod


def main():
    config = configparser.ConfigParser()
    config.read('config_radial_density.ini')

    trajectory_dir = config['DEFAULT']['trajectory_dir']
    output_dir = config['DEFAULT']['output_dir']
    output_filename = config['DEFAULT']['output_filename']
    start_time = int(config['DEFAULT']['start_time'])
    end_time = int(config['DEFAULT']['end_time'])
    interval = int(config['DEFAULT']['interval'])
    water_type = int(config['DEFAULT']['water_type'])
    cutoff_radius = float(config['DEFAULT']['cutoff_radius'])
    bin_width = float(config['DEFAULT']['bin_width'])
    molar_mass_water = float(config['DEFAULT']['molar_mass_water'])
    avogadro_number = float(config['DEFAULT']['avogadro_number'])
    bulk_water_density = float(config['DEFAULT']['bulk_water_density'])

    radial_density_mod.calculate_from_files(
        trajectory_dir,
        output_dir,
        output_filename,
        start_time,
        end_time,
        interval,
        water_type,
        cutoff_radius,
        bin_width,
        molar_mass_water,
        avogadro_number,
        bulk_water_density
    )

if __name__ == '__main__':
    main()
