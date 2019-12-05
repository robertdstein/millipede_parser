import pickle
import os
import numpy as np
from astropy.io import fits
import healpy as hp
import argparse
from parse_archival_scan import get_v0_output_dir, get_v0_output_file

# output_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/output_raw_fits/EHE/"
# cache_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/EHE/"


def parse_archival_txt(candidate, base_output_dir, cache_dir):
    candidate_name = candidate.split(".")[0]
    path = os.path.join(cache_dir, "{0}/{1}.skymap_SpiceMie.txt".format(candidate, candidate_name))
    output_dir = get_v0_output_dir(base_output_dir)

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_name = get_v0_output_file(candidate)

    output_file = os.path.join(output_dir, output_name)

    split = candidate.split("_")
    #
    with open(path, "r") as f:
        data = np.array([float(x) for x in f.readlines()])

    hdu = fits.PrimaryHDU(data=data)
    best_key = list(data).index(min(data))
    best_res = data[best_key]
    hdr = hdu.header
    hdr.set('NSIDE', hp.pixelfunc.npix2nside(len(data)))
    hdr.set('Coord', "ICECUBE_LOCAL")
    # hdr.set('time_mjd', res["time_mjd"])
    # hdr.set('time_utc', res["time_string"])
    # hdr.set("run_id", res["run_id"])
    hdr.set("DATA", "LOGL")
    hdr.set("minpixel", best_key)
    # hdr.set("E_dep", best_res["depositedEnergy"])
    hdr.set("ICEMODEL", "SpiceMie")
    hdr.set("ARCHIVAL", True)
    if "EHE" in split:
        hdr.set("Stream", "EHE")
    elif "HESE" in split:
        hdr.set("Stream", "HESE")
    print("Writing to", output_file)
    hdu.writeto(output_file, overwrite=True)
    return output_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-c", "--cache_dir")
    parser.add_argument("-e", "--event", default=None)
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    else:
        candidates = sorted([y for y in os.listdir(args.cache_dir) if "event" in y])

    for candidate in candidates:
        parse_archival_txt(candidate, args.output_dir, args.cache_dir)