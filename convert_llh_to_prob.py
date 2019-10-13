import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from scipy.stats import chi2, norm
from add_contextual_info import get_v2_output_dir

# ======================================
# The magic Millipede contour numbers
# ======================================

threshold_50 = 22.2
threshold_90 = 64.2

# ======================================


def get_v3_output_dir_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v3_prob_map")


def convert_prob_ts(p):
    return chi2.ppf(p, 2)


def convert_sigma_ts(sig):
    return chi2.ppf(2 * norm.cdf(sig) - 1.0, 2)

expected_lower_50 = convert_prob_ts(0.5)
expected_upper_90 = convert_prob_ts(0.90)
grad_50_90 = (expected_upper_90 - expected_lower_50)/(threshold_90 - threshold_50)

def extrapolate_ts_from_90(p):
    return threshold_90 + (convert_prob_ts(p) - expected_upper_90)/grad_50_90


# Set everything beyond this percentile threshold  to zero
truncation_percentile = 0.99
threshold_max = extrapolate_ts_from_90(truncation_percentile)
contours = [
    (0.0, threshold_50, 0.0, 0.5),
    (threshold_50, threshold_90, 0.5, 0.9),
    (threshold_90, threshold_max, 0.9, truncation_percentile)
]

def apply_mask_rescale(probs, lower_ts, upper_ts, lower_prob, upper_prob):
    mask = np.logical_and(probs > lower_ts,
                          probs < upper_ts
                         )
    probs[mask] -= lower_ts
    expected_lower = convert_prob_ts(lower_prob)
    expected_upper = convert_prob_ts(upper_prob)
    probs[mask] *= (expected_upper - expected_lower)/(upper_ts - lower_ts)
    probs[mask] += expected_lower
    return probs


def convert_to_prob(data):
    """
    Apply four-stage rescale of LLH landscape, to reach Wilk's-like rescaled contours.
    Then converts these to probabilities, using exponential.
    Rescales 0-50%, and 50-90%, based on known contour threshold numbers.
    Extrapolated to 99% contour using 50-90%. Truncates everything beyond 99% to prob of 0.

    :param data: delta-llh landscape
    :return: pixelwise probabilities
    """
    probs = list(np.copy(data))
    probs = np.array(probs) - min(probs)
    probs[np.isnan(probs)] = max(probs)
    for (lower_ts, upper_ts, lower_prob, upper_prob) in contours:
        probs = apply_mask_rescale(probs, lower_ts, upper_ts, lower_prob, upper_prob)
    mask = probs > threshold_max
    probs = np.exp(-probs)
    probs[mask] = 0.0
    probs /= np.sum(probs)
    return probs


def convert_llh_to_prob(candidate, base_output_dir):
    input_dir = get_v2_output_dir(base_output_dir)
    path = os.path.join(input_dir, candidate)
    print(candidate)
    output_dir = get_v3_output_dir_dir(base_output_dir)

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_file = os.path.join(output_dir, candidate)

    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        header["DATA"] = "PROB"
        hdul[0].data = convert_to_prob(data)
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
        candidates = [y for y in os.listdir(get_v2_output_dir(args.output_dir)) if "event" in y]

    for candidate in candidates:
        convert_llh_to_prob(candidate, args.output_dir)