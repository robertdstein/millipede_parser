import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from parse_archival_scan import get_v0_output_dir
from contextual_info import archival_data
import healpy as hp
import logging

def extract_az_zen(nside, index):
    return hp.pix2ang(nside, index)

def add_archival_ehe_info(candidate, data, header):
    # zen, phi = extract_az_zen(header["NSIDE"], header["MINPIXEL"])
    # dec = np.degrees(zen - np.pi/2.)

    stream = header["STREAM"]

    namesplit = [x for x in candidate.split("_")]

    try:
        year = int(namesplit[1])
        num = int(namesplit[4].split(".")[0])

        stream_mask = archival_data.T[4] == stream
        year_mask = archival_data[stream_mask].T[3] == year

        match = archival_data[stream_mask][year_mask][num]

        header.set('time_mjd', match[0].mjd)
        header.set("YEAR", int(match[0].jyear))

    except IndexError:
        pass
    except ValueError:
        pass
    return data, header


def get_v1_output_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v1_with_contextual_info")

def add_contextual_info(candidate, base_output_dir):
    input_dir = get_v0_output_dir(base_output_dir)
    path = os.path.join(input_dir, candidate)
    logging.info(candidate)
    output_dir = get_v1_output_dir(base_output_dir)

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_file = os.path.join(output_dir, candidate)

    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    data, header = add_archival_ehe_info(candidate, data, header)

    hdul = fits.PrimaryHDU(data=data, header=header)

    logging.info("Writing to {0}".format(output_file))
    hdul.writeto(output_file, overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-e", "--event", default=None)
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    else:
        candidates = sorted([y for y in os.listdir(get_v0_output_dir(args.output_dir)) if "event" in y])

    for candidate in candidates:
        add_contextual_info(candidate, args.output_dir)
