import pickle
import os
import numpy as np
import healpy as hp
from astropy.io import fits
import argparse
from parse_archival_scan import get_v0_output_dir, get_v0_output_file

# output_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/output_raw_fits/EHE/"
# cache_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/EHE/"


def parse_archival_fits(candidate, base_output_dir, cache_dir):
    root = os.path.join(cache_dir, candidate)
    if not os.path.isfile(root):
        path = os.path.join(root, "skymap_local.fits")
    else:
        path = root
    print(candidate)
    output_dir = get_v0_output_dir(base_output_dir)

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_name = get_v0_output_file(candidate.split(".")[1] + ".x")
    output_file = os.path.join(output_dir, output_name)

    split = candidate.split("_")
    run_id = split[3][3:]
    event_id = split[4][5:]

    meta_file = os.path.join(root, "properties_local.txt")

    with open(meta_file, "rb") as f:
        for line in f.readlines():
            if "mjd" in line:
                time_mjd = float(line.split(": ")[1])
            elif "(string)" in line:
                time_utc = line.split(": ")[1].split("\n")[0]
            elif "minimum" in line:
                minpixel = int(line.split(" at ")[1])
            elif "energy:" in line:
                best_e = float(line.split("energy: ")[1][:-4])

    with fits.open(path) as hdul:
        data = hdul[1].data["logl"]

    hdu = fits.PrimaryHDU(data=data)
    hdr = hdu.header
    hdr.set('NSIDE', hp.pixelfunc.npix2nside(len(data)))
    hdr.set('time_mjd', time_mjd)
    hdr.set('time_utc', time_utc)
    hdr.set('Coord', "ICECUBE_LOCAL")
    hdr.set("ARCHIVAL", True)
    hdr.set("run_id", run_id)
    hdr.set("event_id", event_id)
    hdr.set("DATA", "LOGL")
    hdr.set("minpixel", minpixel)
    hdr.set("E_dep", best_e)
    hdr.set("Stream", "Diffuse")
    hdr.set("ARCHIVAL", True)
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
        parse_archival_fits(candidate, args.output_dir, args.cache_dir)
