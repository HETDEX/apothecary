#!/usr/bin/env python3
"""
All-sky Mollweide map, EQUATORIAL (RA/Dec), centered on RA = 270 deg:
  - light-gray Milky Way band (galactic |b| < MW_B) underneath
  - red HEALPix pixel points on top (native pixel size)

This version does NOT use hp.mollview / hp.projscatter (which trip a
matplotlib >=3.5 "Normalize instance simultaneously with vmin/vmax" error in
some healpy builds). It renders on matplotlib's own 'mollweide' projection and
uses healpy only for the coordinate transforms:
  - hp.pix2ang(..., lonlat=True)        pixel ID   -> (RA, Dec)
  - hp.Rotator(coord=['C','G'])          equatorial -> galactic (for the band)

Dependencies: healpy, numpy, matplotlib.
"""

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# ------------------------- CONFIG -------------------------
INFILE    = "uniq_hp"    #healpix points
NSIDE     = 32768
NEST      = False        # False = RING
CENTER    = 270.0        # RA at image center
MW_B      = 10.0         # Milky Way band half-width in galactic latitude (deg)
MARKER_S  = 3.0          # arbitrary size for readability, not a real size on sky
BAND_GRAY = "0.82"       # matplotlib gray level (0=black .. 1=white)
OUTFILE   = "healpix_equatorial_c270_mw.png"


def shift(lon_deg):
    """Center on CENTER and flip so RA increases to the LEFT (astro convention).
    Returns radians in [-pi, pi] suitable for a matplotlib 'mollweide' axes."""
    return -np.radians(((lon_deg - CENTER + 180.0) % 360.0) - 180.0)


def main():
    # --- red points: IDs -> (RA, Dec) ---
    npix  = hp.nside2npix(NSIDE)
    raw   = np.loadtxt(INFILE, dtype=np.int64)
    valid = raw[(raw >= 0) & (raw < npix)]                 # drop -999 / out-of-range
    print(f"loaded {raw.size} IDs; kept {valid.size}, dropped {raw.size - valid.size} invalid")
    ra, dec = hp.pix2ang(NSIDE, valid, nest=NEST, lonlat=True)

    # --- Milky Way band: galactic |b| on an equatorial grid ---
    gra  = np.arange(0.0, 360.0, 1.0)
    gdec = np.arange(-90.0, 90.001, 0.5)
    RR, DD = np.meshgrid(gra, gdec)
    _, gb = hp.Rotator(coord=['C', 'G'])(RR.ravel(), DD.ravel(), lonlat=True)
    band_b = np.abs(gb)
    Xg, Yg = shift(RR.ravel()), np.radians(DD.ravel())

    # --- plot on matplotlib's own mollweide axes (no hp.mollview) ---
    fig = plt.figure(figsize=(13, 7))
    ax = fig.add_subplot(111, projection='mollweide')
    ax.grid(True, alpha=0.3, lw=0.5, zorder=0)

    # gray band, UNDER the points. tricontourf triangulates in projected-lon
    # space, so the +/-pi seam is naturally split (band shows on both edges).
    ax.tricontourf(Xg, Yg, band_b, levels=[0.0, MW_B], colors=[BAND_GRAY], zorder=1)

    # red HEALPix points ON TOP
    ax.scatter(shift(ra), np.radians(dec), s=MARKER_S, c='crimson',
               alpha=0.7, edgecolors='none', zorder=5)

    xt = np.radians([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150])
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{int((CENTER - np.degrees(t)) % 360)}\u00b0" for t in xt], fontsize=8)
    ax.set_yticklabels([f"{int(np.degrees(t))}\u00b0" for t in ax.get_yticks()], fontsize=8)
    ax.set_xlabel(f"RA  (centered on {int(CENTER)}\u00b0)")
    ax.set_ylabel("Dec")
    ax.set_title(f"{valid.size:,} HEALPix pixels \u2014 equatorial\n"
                 f"Milky Way band |b|<{int(MW_B)}\u00b0 shaded gray \u2014 NSIDE={NSIDE} "
                 f"{'NEST' if NEST else 'RING'}", pad=20, fontsize=11)

    fig.tight_layout()
    fig.savefig(OUTFILE, dpi=150, bbox_inches='tight')
    print(f"saved {OUTFILE}")


if __name__ == "__main__":
    main()


