# NIRCam Ramp Processing Outline

This file expands the recipe in `README.md` into the processing steps implied by the referenced scripts.

## 0. Data Intake And Checks

1. Download JWST `_uncal.fits` and `_cal.fits` files.
2. Check image footprints and geometry, e.g. with `plot_sci_footprints.py`.
3. Check observation dates manually.

## 1. Align And Stack Calibrated Images

Scripts: `step1_resamp_all_ims.py`, `resamp.py`

1. Find all `j*_cal.fits` files.
2. Read each file's `FILTER` keyword.
3. Group calibrated images by filter.
4. For each filter, submit a Slurm job.
5. Run JWST `TweakRegStep` on that filter's images using Gaia DR3 as the absolute reference catalog.
6. Save per-image `_tweakreg.fits` products.
7. Run JWST `ResampleStep` on the tweaked files.
8. Write stacked/resampled filter products such as `F150W_stacked.fits`.

## 2. Build Initial Candidate Star Catalog

Scripts: `step1--2_find_stars.py`, `step1--2_star_catalog.py`

1. Read a stacked science image.
2. Estimate image scatter using NMAD.
3. Find bright local maxima above `30 * NMAD`.
4. Extract 15 x 15 cutouts around those maxima.
5. Normalize cutouts and median-combine them into an approximate PSF.
6. Write `all_cutouts.fits` and `est_PSF.fits`.
7. Correlate the approximate PSF against each valid image position.
8. Write the correlation image as `my_corr.fits`.
9. Find local maxima in `my_corr.fits` with correlation greater than 0.75.
10. Convert candidate pixel coordinates to RA/Dec using the stacked-image WCS.
11. Write `WD_candidates.ecsv` and `ds9.reg`.

## 3. Reprocess Raw Ramps With Custom Linearity/Flat Handling

Scripts: `step5A_wrap.py`, `step5B_wrap.py`, `step5_nonlin.py`

`step5A_wrap.py` runs one representative `_uncal.fits` for each chip/exposure grouping so reference files are downloaded without many jobs racing each other. `step5B_wrap.py` then batches the remaining `_uncal.fits` files.

The exact JWST pipeline tasks called by `step5_nonlin.py`, in order, are:

1. `GroupScaleStep`
2. `DQInitStep`
3. `SaturationStep`
4. `SuperBiasStep`
5. `RefPixStep`
6. `LinearityStep`
7. `DarkCurrentStep`

After those pipeline calls, `step5_nonlin.py` does custom post-processing:

1. Instantiate `FlatFieldStep` only to retrieve reference files.
2. Retrieve `flat`, `area`, `readnoise`, and `gain` references.
3. Build a correction image from the flat and normalized pixel-area map.
4. Save the calibrated ramp to `_uncallin.fits`.
5. Reopen `_uncallin.fits` in update mode.
6. Divide `SCI` by the correction image.
7. Multiply `SCI` by gain to convert to electrons.
8. Append an `RN` image extension containing read noise in electrons.

`FlatFieldStep.call(...)` is not run directly; the flat-related correction is applied manually.

## 4. Build Empirical PSFs

Scripts: `step4B_PSF_wrap.py`, `PSF_builder_mine.py`

`step4B_PSF_wrap.py` finds `_uncallin.fits` files, groups them by detector/filter, and submits Slurm jobs that call `PSF_builder_mine.py` with the candidate catalog and the relevant images. Note that the current script has the chip/filter loop hard-coded to `["nrcb1_F150W"]`, even though the README says the intended run is one job per chip.

`PSF_builder_mine.py` does the following:

1. Read candidate positions from the input catalog: `x_pix`, `y_pix`, `RA_deg`, and `Dec_deg`.
2. Require each input image to be an `_uncallin.fits` file.
3. Build a design matrix that maps a 2x-oversampled PSF parameter grid onto a 10x-sampled PSF image.
4. For each `_uncallin.fits` image:
   - read the `SCI` ramp cube;
   - mask saturated groups using `GROUPDQ & 2`;
   - assert the expected `TGROUP` value;
   - compute group-time midpoints;
   - read the matching `_tweakreg.fits` WCS and filter name;
   - transform candidate RA/Dec into image pixel coordinates;
   - remove close pairs closer than 15 pixels;
   - drop candidates too close to image edges.
5. Convert the ramp cube into group-to-group differences: `sci[:, 1:] - sci[:, :-1]`.
6. Flatten integrations and groups into a stack of difference images.
7. Build a median difference image.
8. Estimate a 2D background with `Background2D` and `SExtractorBackground`.
9. Subtract that background from both the median difference image and every difference image.
10. Extract 21 x 21 pixel cutout stacks around each surviving star.
11. Iterate twice over rough star selection:
    - median-combine each star's cutout stack;
    - build a rough PSF from the median star images;
    - fit each star's rough flux/amplitude;
    - estimate relative background structure using NMAD;
    - reject stars with unusually structured backgrounds;
    - keep the brightest stars, first up to `2 * KEEP_TOP_N_STARS`, then up to `KEEP_TOP_N_STARS`.
12. Pool cutouts from all input images.
13. Check that enough pooled stars remain.
14. Initialize the PSF model as a Gaussian-like 2x-oversampled grid.
15. Renormalize the PSF so the native-sampled core sum is unity.
16. Fit each star's amplitude, centroid offsets, and sky against the current PSF.
17. Iteratively solve the PSF:
    - refit star amplitudes/centroids/sky;
    - solve the PSF parameters with `miniLM_new`;
    - renormalize the 10x PSF;
    - build model cutouts and residual cutouts;
    - continue for at least 2 iterations and up to 5, stopping when the PSF comparison RMS is small.
18. If verbose mode is enabled, write diagnostic FITS files such as:
    - `the_background.fits`
    - `sci_diffs.fits`
    - `median_diff.fits`
    - `PSF_rough_guess_it=*.fits`
    - `PSF_10x_iter=*.fits`
    - `the_mod_iter=*.fits`
    - `the_resid_iter=*.fits`
    - `the_norm_resid_iter=*.fits`
19. Write the final oversampled PSF as `PSF_10x_<detector>_<filter>.fits`.

## 5. Ramp-Level PSF Photometry

Scripts: `step6_wrap.py`, `step6_do_phot.py`

1. Read the `WD*csv` candidate list.
2. Count usable candidate rows.
3. Find short-wave `_uncallin.fits` files.
4. Split the candidate list across many Slurm jobs.
5. For each target image, load short-wave and paired long-wave ramp cubes.
6. Use the `_tweakreg.fits` WCS files to map candidate RA/Dec into detector coordinates.
7. Mask saturated ramp samples.
8. Form group-to-group difference images.
9. Extract small cutouts around each candidate.
10. Fit the empirical or model PSF to each group-difference cutout.
11. Record per-group short and long photometry, RMS diagnostics, uncertainties, centroids, times, RA, and Dec.
12. Write `photo_subset_*.txt`.
13. Concatenate the subset files into `photo_unflat.txt`.

## 6. Correct Frame-Level Linearity Residuals

Script: `step7_fix_linearity.py`

1. Read `photo_unflat.txt`.
2. For each detector/filter/frame, collect stars with good RMS and sufficient flux.
3. Compare each frame's flux to that star's median flux.
4. Fit spline corrections as a function of source brightness.
5. Write `all_splines.pdf`.
6. Apply corrections to photometry and uncertainties.
7. Write `photo_unflat_linear.txt`.

## 7. Fit And Apply Empirical Flat/Sensitivity Corrections

Scripts: `step7_make_flat.py`, `step7_apply_flat.py`

1. For each filter, read `photo_unflat_linear.txt`.
2. Select bright stars, typically with `min_flux=10000`.
3. Fit detector-relative sensitivity terms and spatial terms.
4. For spline-grid runs such as `S4`, build 2D spline basis terms over detector coordinates.
5. Fit relative sensitivity jointly with per-star true fluxes.
6. Write residual diagnostic plots and JSON fit parameters.
7. Read the JSON sensitivity model.
8. Divide photometry and uncertainties by the fitted relative sensitivity.
9. Write `photo_flattened_linear.txt`.
10. Write `applied_flat.png`.

## 8. Search Light Curves For Candidates

Scripts: `step14_wrap.py`, `step14_find_best_candidates.py`

1. Read `photo_flattened_linear.txt`.
2. Build unique candidate/filter pairs.
3. Submit many Slurm jobs.
4. For each candidate/filter, gather all time-series points.
5. Reject bad points using RMS thresholds.
6. Normalize short and long channels.
7. Fit amplification light-curve models using `amplfication.fits`.
8. Compare the amplification model against a constant-flux model.
9. Keep candidates with sufficient chi-square improvement, consistent short/long amplitudes, and adequate time coverage.
10. Write diagnostic PDFs into `candidate_plots/`.

## 9. Build Median Stellar Flux Catalog

Script: `step8_star_fluxes.py`

1. Read `photo_flattened_linear.txt`.
2. For each star and filter, gather all good measurements with `RMS < 0.2`.
3. Compute median flux per star/filter.
4. Compute uncertainty as median per-point uncertainty divided by `sqrt(N)`.
5. Write `star_fluxes.txt`.

## 10. Match To HST Catalog

Script: `step9_match_to_Sabbi.py`

1. Read `star_fluxes.txt`.
2. Read the PHAT/HST catalog from `../../PHAT_catalogs/v3/merged_catalog.fits`.
3. Cross-match by sky position.
4. Estimate median RA/Dec offset.
5. Apply the offset to the HST catalog.
6. Re-match with a tighter radius.
7. Write `my_with_hst.csv`.

## 11. Prepare Model Atmospheres And AB Conversions

Scripts: `step11_model_atm.py`, `get_AB_abs.py`

1. Read the BOSZ wavelength grid.
2. Read selected BOSZ model atmosphere spectra.
3. Integrate each spectrum through JWST, HST, and Gaia passbands.
4. Compute absolute AB magnitudes for a 1 solar-radius star.
5. Write `model_atmosphere_grid.txt`.
6. Use `get_AB_abs.py` to compute absolute AB magnitudes for Vega, enabling HST Vega magnitudes to be converted to AB.

## 12. Fit Stellar Radii/SEDs

Script: `step10_fit_all_radii_M31.py`

1. Read `my_with_hst.csv`.
2. Read `model_atmosphere_grid.txt`.
3. Convert HST Vega magnitudes to AB.
4. Convert JWST fluxes to AB magnitudes.
5. Fit each star's SED for effective temperature, radius, surface gravity, and extinction.
6. Use the M31 distance modulus.
7. Write chunked outputs like `my_with_hst_fit_XX_of_YY.csv`.
8. Concatenate chunked outputs into `my_with_hst_fit.csv`.

## 13. Bin Stars For Lens-Count Calculations

Script: `step11_get_bins.py`

1. Read `my_with_hst_fit.csv`.
2. Select good SED fits.
3. Bin stars by fitted radius and fractional photometric uncertainty.
4. Convert observation counts into star-hours.
5. For each bin and lens-mass value, submit `step12_get_lens_count.py`.
6. Write `jobs.txt` and `binned_hours.pdf`.
