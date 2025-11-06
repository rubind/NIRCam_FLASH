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
