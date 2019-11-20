import pickle
import os
import numpy as np
from astropy.io import fits
from astropy.time import Time
import argparse
from parse_archival_scan import get_v0_output_dir
from contextual_info import archival_data
import healpy as hp
import logging
import requests

def extract_az_zen(nside, index):
    return hp.pix2ang(nside, index)

def add_archival_ehe_info(candidate, data, header):
    # zen, phi = extract_az_zen(header["NSIDE"], header["MINPIXEL"])
    # dec = np.degrees(zen - np.pi/2.)

    stream = header["STREAM"]

    namesplit = [x for x in candidate.split("_")]

    try:
        year = int(namesplit[1])
        num = int(namesplit[4].split(".")[0])

        stream_mask = archival_data.T[4] == stream
        year_mask = archival_data[stream_mask].T[3] == year

        match = archival_data[stream_mask][year_mask][num]

        header.set('time_mjd', match[0].mjd)
        header.set("YEAR", int(match[0].jyear))

    except IndexError:
        pass
    except ValueError:
        pass
    return data, header

def v1_gcn_url(event_id, run_id):
    return "https://gcn.gsfc.nasa.gov/notices_amon/{0}_{1}.amon".format(int(event_id), int(run_id))

def v2_gcn_url(event_id, run_id):
    return "https://gcn.gsfc.nasa.gov/notices_amon_g_b/{0}_{1}.amon".format(int(run_id), int(event_id))

def retrieve_v2_alert_info(data, header):
    if not header["ARCHIVAL"]:
        url = v2_gcn_url(header["event_id"], header["run_id"])
        page = requests.get(url)
        print("Found GCN: {0}".format(url))

        if not "404 Not Found" in page.text:

            for line in page.text.splitlines():

                row = [x for x in line.split(":")]
                if len(row) > 1:
                    val = [x for x in row[1].split(" ") if x not in ""][0]
                    if row[0] == "ENERGY":
                        header.set("E_TeV", float(val))
                    elif row[0] == "SIGNALNESS":
                        header.set("P_astro", float(val))
                    elif row[0] == "RUN_NUM":
                        header.set("run_id", float(val))
                    elif row[0] == "EVENT_NUM":
                        header.set("event_id", float(val))
                    elif row[0] == "NOTICE_TYPE":
                        if "Gold" in line.split(" "):
                            header.set("Stream", "Gold")
                        elif "Bronze" in line.split(" "):
                            header.set("Stream", "Bronze")
                        else:
                            raise Exception("Stream not found in {0}".format(row))
                    elif row[0] == "FAR":
                        header.set("FAR", float(val))
                    elif row[0] == "DISCOVERY_DATE":
                        disc_date = "20" + line.split(" ")[-2].replace("/", "-")
                    elif row[0] == "DISCOVERY_TIME":
                        disc_time = line.split(" ")[-2][1:-1]

            time_str = "{0} {1} UTC".format(disc_date, disc_time)
            header.set("TIME_UTC", time_str)
            header.set("TIME_MJD", Time("{0}T{1}".format(disc_date, disc_time), scale="utc", format="isot").mjd)

    return data, header

def retrieve_v1_alert_info(data, header):

    if not header["ARCHIVAL"]:
        url = v1_gcn_url(header["event_id"], header["run_id"])
        page = requests.get(url)
        print("Found GCN: {0}".format(url))

        if not "404 Not Found" in page.text:

            for line in page.text.splitlines():

                row = [x for x in line.split(":")]
                if len(row) > 1:
                    val = [x for x in row[1].split(" ") if x not in ""][0]
                    if row[0] == "ENERGY":
                        header.set("E_TeV", float(val))
                    elif row[0] == "SIGNALNESS":
                        header.set("P_astro", float(val))
                    elif row[0] == "RUN_NUM":
                        header.set("run_id", float(val))
                    elif row[0] == "EVENT_NUM":
                        header.set("event_id", float(val))
                    elif row[0] == "NOTICE_TYPE":
                        if "EHE" in line.split(" "):
                            header.set("Stream", "EHE")
                        elif "HESE" in line.split(" "):
                            header.set("Stream", "HESE")
                        else:
                            raise Exception("Stream not found in {0}".format(row))
                    elif row[0] == "FAR":
                        header.set("FAR", float(val))
                    elif row[0] == "DISCOVERY_DATE":
                        disc_date = "20" + line.split(" ")[-2].replace("/", "-")
                    elif row[0] == "DISCOVERY_TIME":
                        disc_time = line.split(" ")[-2][1:-1]

            time_str = "{0} {1} UTC".format(disc_date, disc_time)
            header.set("TIME_UTC", time_str)
            header.set("TIME_MJD", Time("{0}T{1}".format(disc_date, disc_time), scale="utc", format="isot").mjd)

            print(header)



    return data, header


def get_v1_output_dir(base_output_dir):
    return os.path.join(base_output_dir, "fits_v1_with_contextual_info")

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

    data, header = retrieve_v1_alert_info(data, header)
    data, header = retrieve_v2_alert_info(data, header)

    data, header = add_archival_ehe_info(candidate, data, header)

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
        candidates = sorted([y for y in os.listdir(get_v0_output_dir(args.output_dir)) if "event" in y])

    for candidate in candidates:
        add_contextual_info(candidate, args.output_dir)
