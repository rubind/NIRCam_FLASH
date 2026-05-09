Step 0: Download _uncal.fits and _cal.fits files. 

Check geometry: python ~/nic_python/plot_sci_footprints.py *nrcalong_cal.fits --out footprints_MAST_2026-03-31T03_54_35.506Z.png

Check dates, too

python ~/NIRCam_ramp/step1_resamp_all_ims.py (aligns images to Gaia, making _tweakreg.fits files then stacks each band, running more than one of the jobs at the same time may cause reference-file conflicts)

python ~/NIRCam_ramp/step1--2_find_stars.py F150W_stacked_resample.fits (computes approximate PSF, then correlation coefficient with that PSF over the image, no I don't know why the built in function for this has problems)

python ~/NIRCam_ramp/step1--2_star_catalog.py F150W_stacked_resample.fits (writes WD_candidates.ecsv)

python ~/NIRCam_ramp/step5A_wrap.py (runs step5_nonlin.py once for each detector/filter, this will download the reference files without causing conflicts between jobs)

python ~/NIRCam_ramp/step5B_wrap.py

(check if all _uncallin.fits and _tweakreg.fits files got made)

python ~/NIRCam_ramp/step4B_PSF_wrap.py WD_candidates.ecsv (submits 10 jobs, one for each chip)

Modify step6_wrap.sh, which calls this:

python /home/drubin/NIRCam_ramp/step6_wrap.py 0 15 500 0

After running the photometry, concatenate to "photo_unflat.txt." Then:

python ../NIRCam_ramp/step7_fix_linearity.py

for NRCFILT in F150W F277W
do
    python ../NIRCam_ramp/step7_make_flat.py photo_unflat_linear.txt $NRCFILT 10000 0
    python ../NIRCam_ramp/step7_make_flat.py photo_unflat_linear.txt $NRCFILT 10000 S4
done

Now fit light curves to find candidates:

sbatch step14.sh which will do:

cd /home/drubin/supernova_koastore/JWST_programs/MAST_2026-03-31T03_54_35.506Z/JWST

python /home/drubin/NIRCam_ramp/step14_wrap.py 400 photo_flattened_linear.txt

Fitting star radii:

Construct median JWST fluxes and flux uncertainties for all stars:

python step8_star_fluxes.py (reads photo_flattened_linear.txt)

If the model-atmosphere grid needs updating, need to download r=500 grid from https://archive.stsci.edu/hlsp/bosz
