from jwst.resample import ResampleStep
import sys

output_file = sys.argv[1]
these_fls = sys.argv[2:]

ResampleStep.call(these_fls, output_file= output_file)
