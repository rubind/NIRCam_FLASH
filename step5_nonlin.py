from jwst.datamodels import RampModel
from jwst.group_scale import GroupScaleStep
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.refpix import RefPixStep
from jwst.linearity import LinearityStep
import crds
from crds.client import getreferences
import sys

def do_flat():
      refs = get_crds_refs_from_model(ramp)
    flat_path = refs.get("flat")
    area_path = refs.get("area")
    if not flat_path or not os.path.exists(flat_path):
        raise RuntimeError("CRDS could not provide a FLAT reference for this exposure.")
    if not area_path or not os.path.exists(area_path):
        raise RuntimeError("CRDS could not provide an AREA (pixel-area) reference for this exposure.")

    flat = load_flat_image(flat_path)         # dimensionless
    area = load_area_image(area_path)         # steradians per pixel

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
    #   corr = flat * pam_factor
    # We want to divide the image by corr.
    corr = flat_safe * pam_safe

    diffs_corr = diffs / corr  # still DN/s, now flat + PAM corrected




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




newfl = sys.argv[1].replace("_uncal", "_uncallin")
assert newfl != sys.argv[1]

ramp.save(newfl)
