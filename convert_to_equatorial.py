import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from healpy.rotator import Rotator

def switch_ra_azimuth(phi_in, mjd):
    """Givin MJD, transform azimuth->RA **or** RA->azimuth."""
    # Stolen with credit to Mike Richman/csky
    sidereal_day = 0.997269566 # sidereal day = length * solar day
    sidereal_offset = 2.54199002505 # RA = offset + 2pi * (MJD/sidereal_length)%1  - azimuth
    return (sidereal_offset + 2*np.pi*((mjd/sidereal_day)%1) - phi_in) % (2*np.pi)


def rotate_to_equatorial(data, header):
    rot = Rotator(rot=[180. + np.degrees(switch_ra_azimuth(0.0, header["time_mjd"])), 180.0])
    return rot.rotate_map_pixel(data)

def get_input_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v0_raw")


def parse_ehe(candidate, base_output_dir):
    input_dir = get_input_dir(base_output_dir)
    path = os.path.join(input_dir, candidate)
    print(candidate)
    output_dir = os.path.join(base_output_dir, "fits_v1_equatorial")

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_file = os.path.join(output_dir, candidate)

    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        hdul[0].data = rotate_to_equatorial(data, header)
        header["COORD"] = "EQUATORIAL"
        print("Writing to", output_file)
        hdul.writeto(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-e", "--event", default=None)
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    else:
        candidates = [y for y in os.listdir(get_input_dir(args.output_dir)) if "event" in y]

    for candidate in candidates:
        parse_ehe(candidate, args.output_dir)