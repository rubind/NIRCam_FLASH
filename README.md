Step 0: Download _uncal.fits and _cal.fits files.

python ~/NIRCam_ramp/step1_resamp_all_ims.py (aligns images to Gaia, making _tweakreg.fits files then stacks each band, running more than one of the jobs at the same time may cause reference-file conflicts)
