import pickle
import os
import numpy as np
from astropy.io import fits
import argparse

# output_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/output_raw_fits/EHE/"
# cache_dir = "/Users/robertstein/Realtime_Stuff/alert_archive/EHE/"


def parse_ehe(candidate, output_dir, cache_dir):
    path = os.path.join(cache_dir, "{0}/step10_data.pickle".format(candidate))
    print(candidate)
    print("{0}.fits".format(candidate))
    output_file = os.path.join(output_dir, "{0}.fits".format(candidate))

    split = candidate.split("_")

    with open(path, "r") as f:
        x = pickle.load(f)

    data = np.array([(y["logl"]) for y in x[0]["SpiceMie"]], dtype=np.float)
    hdu = fits.PrimaryHDU(data=data)
    best_key = list(data).index(min(data))
    best_res = x[0]["SpiceMie"][best_key]
    hdr = hdu.header
    res = x[0]["SpiceMie"][-1]
    hdr.set('NSIDE', x[1]["SpiceMie"])
    hdr.set('Coord', "ICECUBE_LOCAL")
    hdr.set('time_mjd', res["time_mjd"])
    hdr.set('time_utc', res["time_string"])
    hdr.set("run_id", res["run_id"])
    hdr.set("DATA", "LOGL")
    hdr.set("minpixel", best_key)
    hdr.set("E_dep", best_res["depositedEnergy"])
    hdr.set("ICEMODEL", "SpiceMie")
    hdr.set("ARCHIVAL", True)
    hdr.set("Stream", split[2])
    hdr.set("YEAR", split[1])
    print("Writing to", output_file)
    hdu.writeto(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-c", "--cache_dir")
    parser.add_argument("-e", "--event", default=None)
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    else:
        candidates = [y for y in os.listdir(args.cache_dir) if "candidate" in y]

    for candidate in candidates:
        parse_ehe(candidate, args.output_dir, args.cache_dir)