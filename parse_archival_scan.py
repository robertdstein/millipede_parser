import pickle
import os
import numpy as np
import healpy as hp
from astropy.io import fits
import argparse

# output_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/output_raw_fits/EHE/"
# cache_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/EHE/"

def get_v0_output_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v0_raw")

def get_v0_output_file(candidate):
    cand_name = ".".join(candidate.split(".")[:-1])
    return "{0}.fits".format(cand_name)


def parse_archival_scan(candidate, base_output_dir, cache_dir):
    root = os.path.join(cache_dir, candidate)
    if not os.path.isfile(root):
        path = os.path.join(root, "step10_data.pickle")
    else:
        path = root
    print(candidate)
    output_dir = get_v0_output_dir(base_output_dir)

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_name = get_v0_output_file(candidate)
    output_file = os.path.join(output_dir, output_name)

    split = candidate.split("_")

    with open(path, "r") as f:
        x = pickle.load(f)

    try:
        data = np.array([(y["logl"]) for y in x[0]["SpiceMie"]], dtype=np.float)
        res = x[0]["SpiceMie"][-1]
        best_key = list(data).index(min(data))
        best_res = x[0]["SpiceMie"][best_key]        
        best_e = best_res["depositedEnergy"]
    except KeyError:
        data_8 = np.array([(x[8][key]["llh"]) for key in x[8].keys()], dtype=np.float)
        up_data_64 = hp.pixelfunc.ud_grade(data_8, 64)
        for (key, item) in x[64].items():
            up_data_64[key] = item["llh"]
        up_data_1024 = hp.pixelfunc.ud_grade(up_data_64, 1024)
        for (key, item) in x[1024].items():
            up_data_1024[key] = item["llh"]    
        data = up_data_1024
        res = x[8][0]
        best_key = list(data).index(min(data))
        best_res = x[1024][best_key]
        best_e = best_res['recoLossesInside']

    hdu = fits.PrimaryHDU(data=data)
    hdr = hdu.header
    hdr.set('NSIDE', hp.pixelfunc.npix2nside(len(data)))
    try:
        hdr.set('time_mjd', res["time_mjd"])
        hdr.set('time_utc', res["time_string"])
        hdr.set("run_id", res["run_id"])
        hdr.set('Coord', "ICECUBE_LOCAL")
        hdr.set("ARCHIVAL", True)
    except KeyError:
        hdr.set('Coord', "EQUATORIAL")
        hdr.set("ARCHIVAL", False)
        split_name = candidate.split(".")
        hdr.set("run_id", split_name[0][3:])
        hdr.set("event_id", split_name[1][3:])
    hdr.set("DATA", "LOGL")
    hdr.set("minpixel", best_key)
    hdr.set("E_dep", best_e)
    hdr.set("ICEMODEL", "SpiceMie")
    if "EHE" in candidate:
        hdr.set("Stream", "EHE")
    elif ".i3.bz2_event0000" in candidate:
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
        parse_archival_scan(candidate, args.output_dir, args.cache_dir)
