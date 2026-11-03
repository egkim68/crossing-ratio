# Protocol: a crossing ratio for comparing two groups on a Likert item

Written 2026-11-01, before the confirmatory runs. Two exploratory sets of replications were run while the design was discussed (seed 2026; one ceiling cell and 300 replications of one symmetric cell). They showed that the sign of the latent mean difference is not a usable benchmark under a ceiling, which led to the benchmark below. The confirmatory runs use a new seed.

## Question
A Likert item scored 1 to K and compared between two groups by its mean assumes equal steps between categories. When does the group conclusion depend on that assumption, and does a simple index detect this better than the existing test of equal spacing?

## The crossing ratio
- Exceedance differences: c_k = P(Y_B > k) - P(Y_A > k), k = 1, ..., K - 1. With scores 1 to K, the mean difference D equals the sum of the c_k.
- Any increasing scoring with steps s_1, ..., s_(K-1) gives the mean difference sum(s_k c_k).
- Crossing ratio (CR) = (sum of |c_k| with the sign of D) / (sum of |c_k| with the opposite sign). The sign of D holds for every scoring whose widest step is less than CR times its narrowest. CR = 1 when D = 0.
- CR is infinite when no c_k has the opposite sign (the answer patterns do not cross). Then no scoring of the categories can reverse the result.
- Inference: stratified multinomial bootstrap by group (B = 999). In each resample, CR is computed with the sign of the original D, so values below 1 mean the sign reversed. The percentile lower 95% bound (5th percentile) is reported.
- Flag: a result is "flagged as fragile" when the lower bound is below r0 = 2 (main), with 1.5 and 3 reported as sensitivity values.

## Benchmark (truth)
- Population CR, computed from the true category probabilities in each cell.
- A cell is "fragile" when population CR < 2 and "safe" otherwise (including infinite CR).
- The latent mean difference is not used as truth. Under a ceiling, a group with a higher latent mean can be lower under every scoring of the categories, and no method based on the categories can recover the latent sign. Cells where this happens are counted and reported descriptively.

## Worked example (for the Introduction or Methods)
- Men (n = 10): counts 0, 0, 2, 8, 0 for answers 1 to 5; mean 3.80, SD 0.42.
- Women (n = 20): counts 1, 5, 1, 0, 13; mean 3.95, SD 1.50.
- Welch t = 0.41, df = 24.2, p = 0.68.
- Exceedance differences (women minus men): -0.05, -0.30, -0.15, 0.65. CR = 0.65 / 0.50 = 1.30. The women's higher mean comes only from the top step. If the lower steps were 1.3 times wider than the top step, the means would be equal.

## Comparison methods
- Welch t test on scores 1 to K (the usual practice).
- Wilcoxon-Mann-Whitney test (normal approximation with tie correction).
- Cumulative probit model with a group location effect: likelihood-ratio test of the group effect.
- Equidistance test: likelihood-ratio test of flexible versus equidistant thresholds, both models with the group effect (ordinal::clm, threshold = "equidistant"); 2 df for K = 5 and 4 df for K = 7.

## Simulation (sim_crossing.R, seed 202)
- Latent model: Y* ~ N(mu_g, sigma_g^2), cut at thresholds tau into K categories. Group A: mu = 0, sigma = 1. Group B: mu = delta, sigma = SD ratio.
- K = 5 and 7.
- Threshold patterns (latent scale of group A):
  - K = 5: (1) equal -1.5, -0.5, 0.5, 1.5; (2) narrow middle -1.2, -0.3, 0.3, 1.2; (3) ceiling -2.0, -1.3, -0.7, -0.2; (4) bunched -1.0, -0.2, -0.1, 1.0.
  - K = 7: (1) equal -1.5, -0.9, -0.3, 0.3, 0.9, 1.5; (2) narrow middle -1.6, -0.9, -0.15, 0.15, 0.9, 1.6; (3) ceiling -2.2, -1.7, -1.2, -0.8, -0.4, -0.1; (4) bunched -1.4, -0.7, -0.05, 0.05, 0.7, 1.4.
- delta = 0, 0.2, 0.5; SD ratio = 1, 1.5, 2.
- n per group = 50, 200, 1,000 (equal groups): 2 x 4 x 3 x 3 x 3 = 216 cells.
- Unequal groups: n_A = 100, n_B = 300, K = 5, all patterns, deltas, and SD ratios: 36 cells.
- 1,000 replications per cell (252,000 samples).

## Outcomes
- Type I error (delta = 0, SD ratio 1) of the Welch t, Wilcoxon-Mann-Whitney, and cumulative probit tests.
- Flag rate of CR by cell; in safe cells it is a false flag rate, in fragile cells it is a detection rate.
- Rejection rate of the equidistance test by cell.
- Among significant Welch t results, the share flagged as fragile.
- Coverage of the bootstrap bound: share of samples with lower bound at or below the population CR (finite population CR only).
- Cells where the latent mean difference and the category-level result disagree in sign (descriptive).

## Hypotheses (K = 5, equal groups, unless stated)
- H1. With delta = 0 and SD ratio 1, the Welch t, Wilcoxon-Mann-Whitney, and cumulative probit tests have Type I error rates between 0.025 and 0.075 in every cell.
- H2. In safe cells with delta = 0.5 and n = 1,000, the CR flag rate is at most 0.05 in every cell.
- H3. With SD ratio 1, delta = 0.5, and n = 1,000, the equidistance test rejects in at least 80% of samples under patterns 2 to 4, although these cells are safe.
- H4. In fragile cells with population CR < 1.5 and n = 1,000, the CR flag rate is at least 0.80 in every cell.
- H5. In fragile cells under pattern 1 (equal thresholds) with n = 1,000, the equidistance test rejects in at most 20% of samples.
- H6. The lower bound covers the population CR (lower bound at or below it) in at least 0.925 of samples in every cell with finite population CR and n >= 200.
Decision: supported if the criterion is met in every cell named; one-sided binomial p values against each threshold are reported for the least favorable cell. Other results are descriptive.

## Real data
- GSS 2022 (NORC) and ESS Round 11, downloaded by the author: Likert items with 4 to 7 categories. Groups: gender, age (under 50 versus 50 and over), degree (yes or no).
- Audit: for every item x grouping, Welch t, CR with its bootstrap bound, and the equidistance test. Main summary: the share of significant t results flagged as fragile.

## Code
- R, run by copy and paste; output folders are created automatically.
- sim_crossing.R (simulation), analysis_crossing.R (rates and hypotheses), make_figures.R (PNG figures, separate from the simulation), audit_survey.R (real data).

---

# Study 2: extended confirmatory run

Written 2026-11-03, after Study 1 (seed 202) was completed and before any Study 2 run. Study 1 results are kept as reported, including the three unsupported hypotheses (H2, H5, H6). Study 2 extends the design and tests a new set of hypotheses (H7 to H13). Some of them restate Study 1 hypotheses in a form informed by Study 1 results, and this is stated for each. Seed 404 (Part A), seed 404 with offset 50,000 (Part B), seed 505 (audit).

## Part A: single items (sim2_crossing.R)
- Latent distributions (mean 0, SD 1 in group A): normal; logistic; skewed (standardized chi-square with 4 df). Group B = delta + SD ratio x Z.
- K = 4, 5, 7, each with four threshold patterns (equal, narrow middle, ceiling, bunched); K = 5 and 7 thresholds as in Study 1; K = 4: (1) -1, 0, 1; (2) -0.5, 0, 0.5; (3) -1.6, -1.0, -0.3; (4) -1.0, -0.9, 0.8.
- delta = 0, 0.2, 0.5; SD ratio = 1, 1.25, 1.5, 2.
- Group sizes: 30, 100, 300, 1,000 per group, and an unequal allocation of 100 and 300.
- 3 x 3 x 4 x 3 x 4 x 5 = 2,160 cells, 1,000 replications each (2,160,000 samples); B = 999.
- Methods as in Study 1. The cumulative probit models are fitted to the aggregated table with clm.fit; the null log-likelihood uses its closed form. Results were checked against the formula interface before the run.
- Population CR computed exactly from the category probabilities. Cell types: null (delta = 0 and SD ratio = 1), fragile (population CR < 2), safe (otherwise).

## Part B: summed scales
- One common factor F (group A: N(0, 1); group B: N(delta, SD ratio^2)); m items, each = sqrt(0.6) F + sqrt(0.4) e, cut at the K = 5 thresholds of pattern 1 (equal) or 3 (ceiling); the scale score is the sum of the m items (0 to 4m).
- m = 1, 3, 5; delta = 0.2, 0.5; SD ratio = 1, 1.5, 2; n = 100, 300, 1,000 per group: 108 cells, 1,000 replications each.
- Population CR of the sum computed exactly by numerical integration over F.
- Outcomes: population CR by m (descriptive question Q1: does summing items make a group comparison less dependent on spacing?), CR flag rate, and Welch t rejection.

## Hypotheses for Study 2 (Part A)
- H7 (restates H1 for the wider design). With delta = 0 and SD ratio = 1, the Welch t and Wilcoxon-Mann-Whitney tests have Type I error rates between 0.025 and 0.075 in every cell with n >= 100 per group (all K, distributions, and patterns, including the unequal allocation).
- H8 (replaces H2, informed by Study 1: the lower bound is cautious just above the cutoff). In cells with population CR >= 10 (including infinite) and n = 1,000, the CR flag rate (lower bound < 2) is at most 0.05.
- H9 (extends H4). In cells with population CR < 1.5 and n >= 300, the CR flag rate is at least 0.80.
- H10 (extends H3). Under the normal latent distribution with SD ratio 1, delta = 0.5, and n = 1,000, the equidistance test rejects in at least 80% of samples under patterns 2 to 4, for every K.
- H11 (replaces H5, informed by Study 1). In fragile cells under pattern 1 with the normal latent distribution and n = 1,000, the equidistance test rejects in at most 50% of samples.
- H12 (replaces H6, informed by Study 1: undercoverage occurred only at population CR = 1). The lower bound covers the population CR in at least 0.925 of samples in every cell with 1 < population CR < Inf and n >= 300.
- H13. At n = 1,000, among significant Welch t results, at least 80% are flagged in cells with population CR < 1.5, and at most 20% are flagged in cells with population CR >= 10.
Decision rules as in Study 1: supported if met in every cell named; one-sided binomial p for the least favorable cell. Flag rates at r0 = 1.5 and 3 are reported descriptively.

## Extended real-data audit (audit2_survey.R, seed 505)
- GSS 2022 (21 items) and GSS 2024 (the 18 of these items asked in 2024), groupings gender, age 50 and over, bachelor's degree or higher; unweighted, with weighted CR (wtssps) as a sensitivity check.
- bfi (25 items) and spi (135 items), adults only (age 18 and over); groupings gender, age 30 and over, college degree or higher (bfi education >= 4; spi education >= 6).
- Outcomes: share of significant Welch t results flagged; equidistance rejection rate; agreement of the flag between GSS 2022 and 2024 for the same item and grouping; a list of all significant comparisons whose answer patterns cross.
