import numpy as np
from astropy.time import Time

def parse_archival_alerts():
    archival_data = []
    stream = None
    with open("contextual_info/catalog_of_alerts.txt", "r") as f:
        for line in f.readlines():

            if "EHE" in line:
                stream = "EHE"
            elif "HESE" in line:
                stream = "HESE"

            if line[0] not in ["#", "\n"]:
                vals = [x for x in line.split(" ") if x not in [""]]
                time = Time(float(vals[0]), format="mjd")
                ra = vals[1]
                dec = vals[3]
                archival_data.append((time, float(ra), float(dec), int(time.jyear), stream))

    return np.array(archival_data)

archival_data = parse_archival_alerts()

def switch_ra_azimuth(phi_in, mjd):
    """Givin MJD, transform azimuth->RA **or** RA->azimuth."""
    # Stolen with credit to Mike Richman/csky
    sidereal_day = 0.997269566 # sidereal day = length * solar day
    sidereal_offset = 2.54199002505 # RA = offset + 2pi * (MJD/sidereal_length)%1  - azimuth
    return (sidereal_offset + 2*np.pi*((mjd/sidereal_day)%1) - phi_in) % (2*np.pi)