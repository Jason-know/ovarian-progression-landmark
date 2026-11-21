# Analysis plan recorded before fitting prognostic models

Date 2026-11-20. This retrospective plan was recorded after checking source fields, eligibility, and available follow-up. The study was not prospectively registered.

## Primary question

Among adults with FIGO III–IV serous tubo-ovarian or primary peritoneal carcinoma alive and observed at 12 months after each source's origin, does recording the month of first progression add discrimination for death in the next 12 months beyond knowing whether progression has already occurred? Compare fixed models fitted in TCGA and transported without refitting to AOCS, UKE and ICON7. Use all-cause mortality. The 12-month prediction horizon from the landmark was chosen because the public ICON7 subset has substantial censoring at 24 months after the landmark.

## Model comparisons

- Baseline predictors are age in decades, stage IV, and suboptimal cytoreduction (greater than 1 cm or the source's suboptimal category). The landmark comparator population is identical for all models.
- Status adds progression at or before the landmark. Timing additionally adds months since first progression, set to zero for individuals with no recorded progression by the landmark. The progression coefficient is anchored at six months after progression by centering this term among progressed individuals only.
- Fit logistic models using inverse probability of censoring weights at the prediction horizon. Do not refit the equation, its scaling, or intercept in external studies. Evaluate by cumulative/dynamic AUC, Brier score, observed-to-expected risk, calibration intercept and slope, and exploratory decision curves. The main contrast is timing minus status AUC in each external cohort. Status minus baseline is the principal secondary contrast. Use cohort-specific independent patient bootstrap CIs with fixed model predictions.
- Report stratified five-fold out-of-fold predictions within TCGA. Compare eight fixed algorithms with the same predictors as a secondary analysis without selecting a deployable model. Keep primary predictors fixed. Show penalization paths and shadow-feature screening as descriptive diagnostics.
- Characterize landmark survival by progression timing and fit cohort-specific Cox associations through 24 months following the landmark, adjusted for age, stage and cytoreduction. Adjust external model interpretation for differences in time origin and endpoints.

## Sensitivity analyses

Evaluate landmarks at 6 and 18 months and outcomes 24 months after the 12-month landmark. Exclude individuals progressing exactly in the reported landmark month. Repeat in standard-treatment patients in ICON7 and UKE. Exclude undocumented early TCGA PFI events from fitting. Repeat without age to permit additional external evidence from GSE32062. Repeat without surgical residual status to quantify sensitivity to complete-case selection. Report all planned analyses even when gains are small, estimates are imprecise or models do not transport well.

## Limits of clinical interpretation

This is prognostic modeling among landmark survivors. Risk at the landmark can inform a future evaluation of counseling or study stratification. Cohort-specific definitions of progression and time zero, calendar-era differences, missing residual tumor dimensions and lack of detailed modern maintenance treatment restrict use. No causal mediation analysis or therapeutic recommendations are planned.
