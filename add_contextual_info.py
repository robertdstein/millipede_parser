import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from parse_archival_scan import get_v0_output_dir
import healpy as hp
import logging

def get_v1_output_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v1_with_contextual_info")

def switch_ra_azimuth(phi_in, mjd):
    """Givin MJD, transform azimuth->RA **or** RA->azimuth."""
    # Stolen with credit to Mike Richman/csky
    sidereal_day = 0.997269566 # sidereal day = length * solar day
    sidereal_offset = 2.54199002505 # RA = offset + 2pi * (MJD/sidereal_length)%1  - azimuth
    return (sidereal_offset + 2*np.pi*((mjd/sidereal_day)%1) - phi_in) % (2*np.pi)

def parse_archival_alerts():
    archival_data = []
    with open("contextual_info/catalog_of_alerts.txt", "r") as f:
        for line in f.readlines():
            if line[0] not in ["#", "\n"]:
                vals = [x for x in line.split(" ") if x not in [""]]
                time = vals[0]
                ra = vals[1]
                dec = vals[3]
                archival_data.append((float(time), float(ra), float(dec)))
                # ra_delta = [float(x) / 2.5 for x in vals[2][1:-1].split(",")]
                # dec_delta = [float(x) / 2.5 for x in vals[4][1:-2].split(",")]

    return np.array(archival_data)

archival_data = parse_archival_alerts()

def extract_az_zen(nside, index):
    return hp.pix2ang(nside, index)

def add_archival_info(data, header):
    zen, phi = extract_az_zen(header["NSIDE"], header["MINPIXEL"])
    dec = np.degrees(zen - np.pi/2.)

    delta_dec = abs(np.array(archival_data.T[2] - dec))

    mask = delta_dec <= min(delta_dec) + 0.1

    matches = archival_data[mask]

    for match in matches:

        converted_ra = np.degrees(switch_ra_azimuth(phi, match[0]))

        delta_ra = abs(match[1] - converted_ra)
        if delta_ra < 0.2:
            logging.info("Match found")
            if "TIME_MJD" in header.keys():
                raise Exception("Multiple matches found for {0}".format(dec))
            header.set('time_mjd', match[0])
        else:
            print(converted_ra, match[1])

    if "TIME_MJD" not in header.keys():
        print("Checked")
        print(matches)
        print("DEC", dec)
        print("Name", header)
        raise Exception("No match found")
    return data, header

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

    data, header = add_archival_info(data, header)

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
        candidates = [y for y in os.listdir(get_v0_output_dir(args.output_dir)) if "event" in y]

    for candidate in candidates:
        add_contextual_info(candidate, args.output_dir)