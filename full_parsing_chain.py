import os
import argparse
import numpy as np
from parse_archival_scan import parse_archival_scan, get_v0_output_file, get_v0_output_dir
from parse_archival_txt import parse_archival_txt
from parse_archival_fits import parse_archival_fits
from convert_to_equatorial import convert_to_equatorial
from add_contextual_info import add_contextual_info
from convert_llh_to_prob import convert_llh_to_prob
from create_plot import create_plot

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-c", "--cache_dir")
    parser.add_argument("-e", "--event", default=None)
    parser.add_argument("-d", "--distribution", default="IC160427A", help="IC160427A, IC170922A or diffuse")
    args = parser.parse_args()

    if args.event is not None:
        candidates = [args.event]
    elif args.cache_dir is not None:
        candidates = sorted([y for y in os.listdir(args.cache_dir) if not np.logical_and(
            np.logical_and("event" not in y, "Event" not in y),
            np.logical_and("run" not in y, "Run" not in y),
        )])
    else:
        candidates = sorted([y for y in os.listdir(get_v0_output_dir(args.output_dir)) if "event" in y])

    for candidate in candidates:
        print(candidate)
        try:
            if args.cache_dir is not None:

                try:

                    try:
                        cand_name = parse_archival_scan(candidate, args.output_dir, args.cache_dir)
                    except IOError:
                        cand_name = parse_archival_txt(candidate, args.output_dir, args.cache_dir)

                except IOError:
                    cand_name = parse_archival_fits(candidate, args.output_dir, args.cache_dir)

            else:
                 cand_name = get_v0_output_file(candidate)

            add_contextual_info(cand_name, args.output_dir)
            convert_to_equatorial(cand_name, args.output_dir)
            sys_filename = convert_llh_to_prob(cand_name, args.output_dir, args.distribution)
            create_plot(sys_filename, args.output_dir)

        except KeyError:
            pass

