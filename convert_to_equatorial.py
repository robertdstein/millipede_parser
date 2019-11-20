import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from healpy.rotator import Rotator
from add_contextual_info import get_v1_output_dir
from contextual_info import switch_ra_azimuth


def rotate_to_equatorial(data, header):
    rot = Rotator(rot=[180. + np.degrees(switch_ra_azimuth(0.0, header["time_mjd"])), 180.0])
    return rot.rotate_map_pixel(data)

def get_v2_output_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v2_equatorial")


def convert_to_equatorial(candidate, base_output_dir):
    input_dir = get_v1_output_dir(base_output_dir)
    path = os.path.join(input_dir, candidate)
    print(candidate)
    output_dir = get_v2_output_dir(base_output_dir)

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_file = os.path.join(output_dir, candidate)

    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        if header["COORD"] == "ICECUBE_LOCAL":
            hdul[0].data = rotate_to_equatorial(data, header)
            header["COORD"] = "EQUATORIAL"
        if header["COORD"] == "ICECUBE_INV":
            hdul[0].data = Rotator(rot=[0., 180.0])
            header["COORD"] = "EQUATORIAL"
        print("Writing to", output_file)
        hdul.writeto(output_file, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-e", "--event", default=None)
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    else:
        candidates = sorted([y for y in os.listdir(get_v1_output_dir(args.output_dir)) if "event" in y])

    for candidate in candidates:
        convert_to_equatorial(candidate, args.output_dir)