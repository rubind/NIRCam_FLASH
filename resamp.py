from jwst.resample import ResampleStep
from jwst.tweakreg import TweakRegStep
import sys

output_file = sys.argv[1]
these_fls = sys.argv[2:]



TweakRegStep.call(
    these_fls,
    # Use an absolute reference catalog:
    abs_refcat='GAIADR3',      # JWST tweakreg’s Gaia DR3 keyword

    # Good starting values for a dense NIRCam stellar field:
    expand_refcat=True,
    fitgeometry='rshift',      # rotation + shift
    searchrad=200.,            # arcsec; for relative matches (image–image)
    use2dhist=True,
    separation=0.5,            # min separation (arcsec) in your catalogs
    tolerance=0.3,             # matching tolerance (arcsec), tweak as needed

    # Absolute (Gaia) matching knobs (optional; defaults are okay-ish):
    abs_searchrad=6.0,         # arcsec
    abs_use2dhist=True,
    abs_fitgeometry='rshift',
    save_results=True
)


tweak_fls = [fl.replace("_cal.fits", "_tweakreg.fits") for fl in these_fls]


ResampleStep.call(tweak_fls, output_file= output_file)
