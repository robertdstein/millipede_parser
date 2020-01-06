import pickle
import os
import numpy as np
from astropy.io import fits
import argparse
from scipy.stats import chi2, norm
from convert_llh_to_prob import get_v3_output_dir, get_systematics_filename
from matplotlib import cm
import healpy as hp
import matplotlib.pyplot as plt

fermi_catalogue = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gll_psc_v19.fit")

def get_plot_output_dir_dir(base_output_dir):
    return os.path.join(base_output_dir, "plots")

def find_pixel_threshold(probs, threshold):
    ranked_pixels = np.sort(probs)[::-1]
    int_sum = 0.0
    pixel_threshold = 0.0

    for i, prob in enumerate(ranked_pixels):
        int_sum += prob
        if int_sum > threshold:
            pixel_threshold = prob
            break
    return pixel_threshold

def get_best_pixel(probs, nside):
    max_index = list(probs).index(max(probs))
    dec, ra = hp.pixelfunc.pix2ang(nside, max_index)
    ra = np.degrees(ra)
    dec = np.degrees(np.pi/2. - dec)
    return ra, dec

# Use Green's theorem to compute the area
# enclosed by the given contour.
def area(vs):
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a


def create_plot(candidate, base_output_dir):
    input_dir = get_v3_output_dir(base_output_dir)
    path = os.path.join(input_dir, candidate)
    
    output_dir = get_plot_output_dir_dir(base_output_dir)
    
    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    output_file = os.path.join(output_dir, "{0}.pdf".format(os.path.splitext(candidate)[0]))
    
    with fits.open(path) as hdul:
        probs = hdul[0].data
        header = hdul[0].header
    
    nside = hp.pixelfunc.npix2nside(len(probs))
    
    threshold_90 = find_pixel_threshold(probs, 0.9)
    mask = probs > (threshold_90 / 5.)
    
    pos = np.array([hp.pixelfunc.pix2ang(nside, i, lonlat=True)
                    for i in np.array(range(len(probs)))[mask]]).T

    min_ra = min(pos[0])
    max_ra = max(pos[0])
    min_dec = min(pos[1])
    max_dec = max(pos[1])

    ra, dec = get_best_pixel(probs, nside)

    wrap_around = False

    if not np.logical_and(ra > min_ra, max_ra > ra):
        min_ra += 360
        wrap_around = True

    ra_delta = abs(max_ra - min_ra)

    ras = np.linspace(min([min_ra, max_ra]), min([min_ra, max_ra]) + ra_delta, 101)
    decs = np.linspace(min_dec, max_dec, 101)
    
    log_p = np.log(probs)

    ps = []

    for d in decs:
        row = []
        for r in ras:
            row.append(np.exp(hp.pixelfunc.get_interp_val(log_p, r, d, lonlat=True)))
        #             row.append(hp.pixelfunc.get_interp_val(probs, r, d, lonlat=True))
        ps.append(row)

    ps = np.array(ps)

    threshold_50 = find_pixel_threshold(probs, 0.5)
    threshold_90 = find_pixel_threshold(probs, 0.9)
    threshold_99 = find_pixel_threshold(probs, 0.99)

    levels = [threshold_90, threshold_50]

    y_inches = 3.0
    x_inches = 6.

    fig = plt.figure(figsize=[x_inches, x_inches])

    ax = plt.subplot(111)
    image = ax.pcolormesh(ras, decs, ps, vmin=0, vmax=max(probs), cmap="inferno")
    plt.scatter(ra, dec, marker="*", color="black")
    CS = ax.contour(ras, decs, ps, colors="white", levels=levels)
    CS_99 = ax.contour(ras, decs, ps, colors="white", alpha=0., levels=[threshold_99, threshold_90, threshold_50])

    all_decs = []
    all_ras = []
    for segs in CS_99.allsegs:
        for k in segs:
            all_decs += list(k.T[1])
            all_ras += list(k.T[0])

    mid_y = np.mean([min(all_decs), max(all_decs)])
    mid_x = np.mean([min(all_ras), max(all_ras)])

    width = max([(max(all_decs) - min(all_decs)), (max(all_ras) - min(all_ras))]) / 2.

    lower_x = mid_x - width
    upper_x = mid_x + width
    lower_y = mid_y - width * abs(np.cos(np.radians(mid_y)))
    upper_y = mid_y + width * abs(np.cos(np.radians(mid_y)))

    ax.set_ylim(lower_y, upper_y)
    ax.set_xlim(lower_x, upper_x)

    contour_labels = ["90%", "50%"]

    npix = (4. * width ** 2) / hp.nside2pixarea(nside, degrees=True)

    for i in range(len(levels)):
        a = 0.
        for m in CS.collections[i].get_paths():
            a += abs(area(m.vertices))
        a *= abs(np.cos(np.radians(mid_y)))
        CS.collections[i].set_label('{0} - area: {1:.2f} sqdeg'.format(contour_labels[i], a))

    ax.clabel(CS, inline=1, fontsize=10, fmt={threshold_99: "99%", threshold_90: "90%", threshold_50: "50%"})
    ax.scatter(ra, dec, marker="*", color="black", s=100)

    with fits.open(fermi_catalogue) as hdul:
        fgl = hdul[1].data

    if wrap_around:
        fgl['RAJ2000'][fgl['RAJ2000'] < lower_x] += 360.

    mask = np.logical_and(
        np.logical_and(
            fgl['RAJ2000'] > lower_x,
            fgl['RAJ2000'] < upper_x
        ),
        np.logical_and(
            fgl['DEJ2000'] > lower_y,
            fgl['DEJ2000'] < upper_y
        )
    )

    if np.sum(mask) > 0:

        ax.scatter(
            fgl['RAJ2000'][mask], fgl['DEJ2000'][mask], s=50, marker='o', color='orange'
        )

        ax.scatter(
            fgl['RAJ2000'][mask], fgl['DEJ2000'][mask], s=10, marker='o', color='navy'
        )

        if np.sum(mask) < 30:
            for i in range(np.sum(mask)):
                ax.text(
                    fgl['RAJ2000'][mask][i] - 0.05 * width, fgl['DEJ2000'][mask][i],
                    fgl['Source_Name'][mask][i], color='white', bbox=dict(facecolor='black', alpha=0.5),
                    fontsize=10
                )

    ax.invert_xaxis()
    ax.grid(which="both", linestyle=":")
    ax.legend()
    ax.set_facecolor('black')
    ax.set_xlabel("Right Ascension (deg)")
    ax.set_ylabel("Declination (deg)")
    plt.savefig(output_file)
    plt.close()
    # fig = hp.gnomview(probs, rot=[mid_x, mid_y], xsize=npix)
    # fig = hp.gnomview(np.log(probs), rot=[mid_x, mid_y], xsize=npix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-e", "--event", default=None)
    parser.add_argument("-d", "--distribution", default="IC160427A", help="IC160427A, IC170922A or diffuse")
    args = parser.parse_args()
    
    if args.event is not None:
        candidates = [get_systematics_filename(args.event, args.distribution)]
    
    else:
        candidates = sorted([y for y in os.listdir(get_v3_output_dir(args.output_dir)) if "event" in y])

    for candidate in candidates:
        create_plot(candidate, args.output_dir)
