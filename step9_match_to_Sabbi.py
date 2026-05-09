import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

def xmatch_nearest_within_radius(
    cat1: pd.DataFrame,
    cat2: pd.DataFrame,
    ra1="RA", dec1="Dec",
    ra2="RA", dec2="Dec",
    radius_arcsec=0.2,
    suffix2="_hst",
):
    """
    For each source in cat1, find nearest source in cat2 within radius.
    Returns cat1 with matched cat2 columns appended (suffix2), plus match distance.
    Unmatched rows get NaNs in appended columns.
    """
    c1 = SkyCoord(ra=cat1[ra1].to_numpy()*u.deg, dec=cat1[dec1].to_numpy()*u.deg, frame="icrs")
    c2 = SkyCoord(ra=cat2[ra2].to_numpy()*u.deg, dec=cat2[dec2].to_numpy()*u.deg, frame="icrs")

    idx, sep2d, _ = c1.match_to_catalog_sky(c2)  # nearest neighbor in cat2 for each cat1 row

    matched = sep2d <= (radius_arcsec * u.arcsec)

    out = cat1.copy()
    out["match_dist_arcsec"] = sep2d.to(u.arcsec).value
    out["matched"] = matched

    # Append cat2 columns (with suffix) for matched rows; NaN for unmatched
    cat2_renamed = cat2.copy()
    cat2_renamed = cat2_renamed.add_suffix(suffix2)

    # Prepare an empty frame with same length as cat1 to merge in
    add = pd.DataFrame(index=np.arange(len(cat1)))
    for col in cat2_renamed.columns:
        add[col] = np.nan

    # Fill matched rows
    add.loc[matched, :] = cat2_renamed.iloc[idx[matched]].to_numpy()

    out = pd.concat([out.reset_index(drop=True), add.reset_index(drop=True)], axis=1)

    RA_offset = np.nanmedian(  (out["RA"] - out["RA_hst"])  )
    Dec_offset = np.nanmedian(  (out["Dec"] - out["Dec_hst"])  )

    print("RA_offset", RA_offset*3600.)
    print("Dec_offset", Dec_offset*3600.)
    return out, RA_offset, Dec_offset

# --- Example usage ---
my = pd.read_csv("star_fluxes.txt", sep=r"\s+")
print(my)
#hst = pd.read_csv("hlsp_http_hst_acs-wfc3_tarantula_multi_v2.0_cat.txt", sep=r"\s+")
hst_tab = Table.read("some_catalog.fits", hdu=1)
hst = hst_tab.to_pandas()



print(hst)
out, RA_offset, Dec_offset = xmatch_nearest_within_radius(my, hst, ra1="RA", dec1="Dec", ra2="ra", dec2="dec", radius_arcsec=0.25)

hst["ra"] += RA_offset
hst["dec"] += Dec_offset

out, NA, NA = xmatch_nearest_within_radius(my, hst, ra1="RA", dec1="Dec", ra2="RA", dec2="Dec", radius_arcsec=0.1)
out.to_csv("my_with_hst.csv", index=False)

