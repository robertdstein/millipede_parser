import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from convert_to_equatorial import get_v1_output_dir

def get_v2_output_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v2_with_contextual_info")


def add_contextual_info(candidate, base_output_dir):
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

        # Insert magic here

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
        candidates = [y for y in os.listdir(get_v1_output_dir(args.output_dir)) if "event" in y]

    for candidate in candidates:
        add_contextual_info(candidate, args.output_dir)