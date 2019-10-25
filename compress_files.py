import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from scipy import sparse
from convert_llh_to_prob import get_v3_output_dir

def get_compressed_output_dir(base_output_dir):
    return os.path.join(base_output_dir, "compressed_files")

def compress_files(candidate, base_output_dir):
    input_dir = get_v3_output_dir(base_output_dir)
    path = os.path.join(input_dir, candidate)
    output_dir = get_compressed_output_dir(base_output_dir)

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    # Here the data is saved as a sparse array, and the header info is saved to the pickle file

    output_npz_file = os.path.join(output_dir, "{0}.npz".format((os.path.splitext(candidate)[0])))
    output_pkl_file = os.path.join(output_dir, "{0}.pkl".format((os.path.splitext(candidate)[0])))

    with fits.open(path) as hdul:
        probs = hdul[0].data
        header = hdul[0].header

    mask = probs > 0.

    sparse.save_npz(output_npz_file, sparse.csr_matrix(mask, dtype=np.bool))

    res_dict = dict()

    for x in header:
        res_dict[x] = header[x]

    with open(output_pkl_file, "wb") as f:
        res_dict["prob"] = probs[mask]
        res_dict["output_path"] = output_npz_file
        pickle.dump(res_dict, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-e", "--event", default=None)
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    else:
        candidates = [y for y in os.listdir(get_v3_output_dir(args.output_dir)) if "event" in y]

    for candidate in candidates:
        compress_files(candidate, args.output_dir)