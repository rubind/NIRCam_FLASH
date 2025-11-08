from jwst.datamodels import RampModel
from jwst.datamodels import open as dm_open
from jwst.group_scale import GroupScaleStep
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.dark_current import DarkCurrentStep
from jwst.refpix import RefPixStep
from jwst.linearity import LinearityStep
from jwst.datamodels import ReferenceFileModel
from jwst.flatfield import FlatFieldStep
import numpy as np
from astropy.io import fits

#import crds
#from crds.client import getreferences
#from crds.core import rmap
import sys

def load_flat_image(flat_path):
    """Open the pipeline FLAT reference and return the combined flat array (float32)."""
    # FLAT refs are delivered as an ImageModel-like file with .data (combined S×D×F flat)
    with dm_open(flat_path) as fm:
        flat = np.array(fm.data, dtype=np.float32)
    return flat

def load_area_image(area_path):
    """Open the AREA reference and return the pixel-area (steradians) array (float32)."""
    with dm_open(area_path) as am:
        area = np.array(am.data, dtype=np.float32)  # sr/pixel
    return area


def do_flat(ramp):
    #refs = get_crds_refs_from_model(ramp)
    #flat_path = refs.get("flat")
    #area_path = refs.get("area")

    #flat_ref = ramp.get_reference_file('flat')
    #area_ref = ramp.get_reference_file('area')

    
    #if not flat_path or not os.path.exists(flat_path):
    #    raise RuntimeError("CRDS could not provide a FLAT reference for this exposure.")
    #if not area_path or not os.path.exists(area_path):
    #    raise RuntimeError("CRDS could not provide an AREA (pixel-area) reference for this exposure.")
    step = FlatFieldStep()  # any Step works; use the one relevant to the ref
    flat_ref = step.get_reference_file(ramp, 'flat')   # path or 'N/A'
    area_ref = step.get_reference_file(ramp, 'area')   # path or 'N/A'
    readnoise_ref = step.get_reference_file(ramp, 'readnoise')
    gain_ref = step.get_reference_file(ramp, 'gain')   # 'gain' is the REFTYPE

    if flat_ref in (None, 'N/A') or area_ref in (None, 'N/A'):
        raise RuntimeError(f"Missing refs: flat={flat_ref}, area={area_ref}")


    read_noise = load_flat_image(readnoise_ref)
    gain_eminus_per_ADU = load_flat_image(gain_ref)
    median_eminus_per_ADU = np.nanmedian(gain_eminus_per_ADU)

    assert (median_eminus_per_ADU > 1)*(median_eminus_per_ADU < 4)
    
    flat = load_flat_image(flat_ref)         # dimensionless
    area = load_area_image(area_ref)         # steradians per pixel
    data = ramp.data.astype(np.float32)
    print("data", data.shape)
    nints, ng, ny, nx = data.shape

    if flat.shape != (ny, nx):
        raise RuntimeError(f"Flat shape {flat.shape} != image {(ny, nx)}")
    if area.shape != (ny, nx):
        raise RuntimeError(f"AREA shape {area.shape} != image {(ny, nx)}")

    # Build the "photometry-ready" correction:
    #  - Apply flat (divide by flat)
    #  - Apply pixel-area correction so point-source photometry is position-independent.
    #
    # We use an area *normalization* so the correction is dimensionless and keeps units in DN/s.
    # Define pam_factor = AREA / median(AREA). Then divide image by pam_factor.
    # => effectively multiply by (median(AREA)/AREA).
    pam_factor = area / np.nanmedian(area)
    # Avoid divide-by-zero / invalids
    eps = 1e-8
    good = np.isfinite(flat) & (flat > eps) & np.isfinite(pam_factor) & (pam_factor > eps)
    flat_safe = flat.copy()
    pam_safe = pam_factor.copy()
    flat_safe[~good] = np.nan
    pam_safe[~good] = np.nan

    # Combined multiplicative correction map:
    corr = flat_safe / pam_factor
    # We want to divide the image by corr.
    #corr = flat_safe #* pam_safe No PAM, as flats are supposed to be uniform detector illumination?
    return corr, read_noise, gain_eminus_per_ADU


# Load your uncal ramp
assert sys.argv[1].count("_uncal.fits") == 1
ramp = RampModel(sys.argv[1])



# Run the subset of steps needed so LINEARITY is applied correctly
ramp = GroupScaleStep.call(ramp)   # e-/DN scaling if needed
ramp = DQInitStep.call(ramp)
ramp = SaturationStep.call(ramp)   # populates SAT DQ flags using reference full-well
ramp = SuperBiasStep.call(ramp)    # subtract superbias
ramp = RefPixStep.call(ramp)       # reference pixel correction
#ramp = FirstFrameStep.call(ramp)      # correct reset anomaly on first group
#ramp = LastFrameStep.call(ramp)       # correct last-group anomaly
ramp = LinearityStep.call(ramp)    # <-- applies per-pixel polynomial linearization to groups
ramp = DarkCurrentStep.call(ramp)     # subtract dark current + amp glow using ref files


# Now ramp.data contains linearized groups (shape: nints, ngroups, ny, nx)

corr, read_noise, gain_eminus_per_ADU = do_flat(ramp)



newfl = sys.argv[1].replace("_uncal", "_uncallin")
assert newfl != sys.argv[1]

ramp.save(newfl)

f = fits.open(newfl, 'update')

print("Before flat", np.sqrt(np.nanmean(np.square(f["SCI"].data))))
f["SCI"].data /= corr
f["SCI"].data *= gain_eminus_per_ADU

f.append(fits.ImageHDU(data=read_noise, name="RN"))
print("After flat", np.sqrt(np.nanmean(np.square(f["SCI"].data))))


f.flush()
f.close()
