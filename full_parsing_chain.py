import os
import argparse
from parse_archival_scan import parse_archival_scan
from parse_archival_txt import parse_archival_txt
from convert_to_equatorial import convert_to_equatorial
from add_contextual_info import add_contextual_info
from convert_llh_to_prob import convert_llh_to_prob
from create_plot import create_plot

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
        try:
            cand_name = parse_archival_scan(candidate, args.output_dir, args.cache_dir)
        except IOError:
            cand_name = parse_archival_txt(candidate, args.output_dir, args.cache_dir)
        add_contextual_info(cand_name, args.output_dir)
        convert_to_equatorial(cand_name, args.output_dir)
        convert_llh_to_prob(cand_name, args.output_dir)
        create_plot(cand_name, args.output_dir)
