import numpy as np
from astropy.time import Time

# Gathered for i3live

ic_season_start_dates = {
    2010: Time("2010-05-31T07:08:28", format="isot", scale="utc"),
    2011: Time("2011-05-13T09:59:43", format="isot", scale="utc"),
    2012: Time("2012-04-26T10:09:18", format="isot", scale="utc"),
    2013: Time("2013-04-18T11:46:20.3992707716", format="isot", scale="utc"),
    2014: Time("2014-04-10T09:23:35", format="isot", scale="utc"),
    2015: Time("2015-04-24T01:42:40", format="isot", scale="utc"),
    2016: Time("2016-05-20T20:50:02.7960466256", format="isot", scale="utc"),
    2017: Time("2017-05-18T04:04:40.3975572356", format="isot", scale="utc"),
    2018: Time("2018-07-10T17:52:03", format="isot", scale="utc")
}


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

                year = None

                for y in sorted(ic_season_start_dates.keys()):
                    if year is None:
                        if np.logical_and(
                                time > ic_season_start_dates[y],
                                time < ic_season_start_dates[y+1]):
                            year = y

                archival_data.append((time, float(ra), float(dec), year, stream))

    return np.array(archival_data)

archival_data = parse_archival_alerts()

def switch_ra_azimuth(phi_in, mjd):
    """Givin MJD, transform azimuth->RA **or** RA->azimuth."""
    # Stolen with credit to Mike Richman/csky
    sidereal_day = 0.997269566 # sidereal day = length * solar day
    sidereal_offset = 2.54199002505 # RA = offset + 2pi * (MJD/sidereal_length)%1  - azimuth
    return (sidereal_offset + 2*np.pi*((mjd/sidereal_day)%1) - phi_in) % (2*np.pi)