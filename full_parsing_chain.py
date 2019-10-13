import os
import argparse
from parse_archival_scan import parse_archival_scan
from convert_to_equatorial import convert_to_equatorial
from add_contextual_info import add_contextual_info

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-c", "--cache_dir")
    parser.add_argument("-e", "--event", default=None)
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    else:
        candidates = [y for y in os.listdir(args.cache_dir) if "event" in y]

    for candidate in candidates:
        cand_name = parse_archival_scan(candidate, args.output_dir, args.cache_dir)
        convert_to_equatorial(cand_name, args.output_dir)
        add_contextual_info(cand_name, args.output_dir)
