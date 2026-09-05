# Manuscript framing and replication plan

## Core claim

> **The genome predicts susceptibility; a sufficiently informative clinical history predicts disease. Once the latter becomes decisive, the former adds little.**

The manuscript should be framed around a general statement about **regimes of predictive information**, not around ZeBRA outperforming genomics per se.

For a non-Mendelian disease, germline genomic information provides a static, finite likelihood contribution. Longitudinal clinical history provides an accumulating information channel. When a near-term clinical predictor derived from that history reaches a sufficiently high-discrimination / high-likelihood-ratio regime, the remaining incremental value of genomic information becomes localized to the predictor's decision boundary. In the extreme clinical-evidence tail, bounded genomic evidence cannot materially change the optimal decision.

The intended claim is therefore stronger and more specific than the generic statement that EHR and genomic predictors are complementary:

- genomics can be important for susceptibility and long-horizon prediction;
- genomics can remain useful where clinical evidence is ambiguous;
- but **for sufficiently accurate near-term clinical prediction, genomic information should add little globally and should matter primarily near the operating boundary**;
- at extreme specificity, the clinical likelihood-ratio tail should dominate bounded genomic evidence.

## Contrast with Detrois et al.

A key recent comparison paper is:

K. E. Detrois, T. Hartonen, M. Teder-Laving, et al., **"Cross-biobank generalizability and accuracy of electronic health record-based predictors compared to polygenic scores,"** *Nature Genetics*, vol. 57, pp. 2136–2145, 2025. DOI: **10.1038/s41588-025-02298-9**.

Detrois et al. trained EHR-derived phenotype risk scores for 13 diseases across FinnGen, UK Biobank, and the Estonian Biobank and reported that EHR-derived scores and PGS captured largely independent information and often provided additive predictive benefit when combined.

Our manuscript should not dispute the empirical result in the regime they studied. The stronger point is that **"EHR + genetics are complementary" is not a universal property of the two modalities**. Complementarity should depend on prediction horizon, clinical predictor strength, and position in clinical likelihood-ratio space.

The distinction to emphasize:

**Detrois et al.:** Are EHR-derived risk scores and PGS complementary over relatively long-horizon disease-onset prediction?

**Our question:** Once a near-term longitudinal clinical predictor becomes highly informative, where can bounded genomic evidence still improve prediction?

**Our answer:** Primarily near the decision boundary; not throughout risk space, and progressively less in the extreme clinical-evidence tail.

This should be expressed carefully. We should not claim that genomics is biologically unimportant. The claim is about **incremental predictive utility conditional on a sufficiently strong clinical predictor**.

## Current empirical anchor: fibrotic ILD / IPF

The Colorado ILD/IPF analysis already provides the first empirical example.

Key observations to retain:

- ZeBRA is the substantially stronger single predictor for the near-term fibrotic ILD/FILA endpoint.
- A broad selected genomic feature panel does not materially improve global ZeBRA AUC.
- MUC5B (`rs35705950`) has a strong and biologically established disease association.
- MUC5B information becomes useful in restricted ZeBRA score regions rather than across the entire risk distribution.
- The fixed regional MUC5B rule produced improved sensitivity near a selected operating boundary at approximately unchanged FPR, despite no global AUC improvement.
- This is qualitatively consistent with the likelihood-ratio boundary-localization theorem.

The paper should explicitly distinguish:

1. **global discrimination** — little or no improvement from adding genomics to a strong near-term ZeBRA predictor; and
2. **local decision utility** — measurable genomic value in a restricted region near the clinical decision boundary.

## Planned disease replication: ADRD

The next high-priority analysis is to reproduce the complete result in Alzheimer's disease / ADRD, where there is also a strong and well-characterized genomic susceptibility architecture.

### ADRD analysis checklist

- [ ] Define the same type of near-term ADRD prediction endpoint used by the published ZeBRA ADRD model.
- [ ] Reconstruct the patient-level clinical/genomic cohort using the same cohort discipline used in the corrected ILD analysis.
- [ ] Include **APOE** genotype explicitly.
- [ ] Include an established ADRD/AD polygenic risk score or the strongest defensible available genomic predictor, rather than relying only on selected individual SNPs.
- [ ] Evaluate genomic-only, ZeBRA-only, and combined predictors under identical held-out splits.
- [ ] Report repeated-split AUC distributions and paired changes.
- [ ] Evaluate calibration and proper scoring rules in addition to AUC.
- [ ] Estimate local genomic utility as a function of ZeBRA percentile / clinical likelihood-ratio region.
- [ ] Identify whether APOE/PRS information is enriched among ZeBRA false negatives or other boundary-adjacent cases.
- [ ] Construct prespecified or training-derived genomic rescue / switching rules and evaluate them on held-out data.
- [ ] Compare sensitivity, specificity, PPV, NPV, LR+, and LR- at matched FPR operating points.
- [ ] Test whether genomic gain decreases in the extreme high-specificity clinical tail.
- [ ] Determine whether the same qualitative pattern seen for MUC5B/IPF is reproduced: **minimal global gain but localized boundary gain**.

The strongest result would be the same qualitative geometry in two very different diseases with prominent genetic drivers:

- **MUC5B / fibrotic ILD-IPF**
- **APOE + ADRD PRS / Alzheimer's disease-ADRD**

This would make it difficult to dismiss the ILD finding as disease-specific.

## Independent cohort replication: All of Us

A second major step is independent replication in the **All of Us Research Program**.

Ideal design:

| Disease | Colorado genomic-clinical cohort | All of Us |
| --- | --- | --- |
| Fibrotic ILD / IPF | current analysis | planned replication |
| ADRD | planned analysis | planned replication |

### All of Us checklist

- [ ] Confirm availability and usable representation of the relevant germline variants / genomic data for MUC5B and ADRD-associated loci.
- [ ] Construct disease-specific genomic predictors using only information available within the cohort and/or externally defined PRS weights.
- [ ] Reconstruct ZeBRA-compatible longitudinal EHR histories with matched prediction horizons.
- [ ] Reproduce the same target definitions as closely as the data permit.
- [ ] Use a strictly external replication design: freeze the analysis logic and boundary tests before inspecting All of Us results.
- [ ] Replicate global ZeBRA vs genomic vs combined discrimination.
- [ ] Replicate local/boundary-conditioned genomic utility.
- [ ] Replicate matched-FPR operating-point analyses.
- [ ] Replicate the extreme-specificity test predicted by the theorem.

If both diseases reproduce in both datasets, the empirical design becomes a 2 x 2 validation of the general principle rather than a single-disease observation.

## Horizon-dependent analysis

A particularly important extension is to vary prediction horizon.

For each disease, where feasible evaluate horizons such as:

- [ ] 1 year
- [ ] 2 years
- [ ] 5 years
- [ ] 10 years

The expected qualitative result is that genomic information should be relatively more useful at longer horizons, while the clinical-history channel should dominate increasingly as the event becomes near-term and the clinical predictor becomes more discriminative.

This provides a direct bridge between the Detrois et al. regime and the ZeBRA regime:

- long-horizon, weaker clinical evidence -> more room for genomic complementarity;
- near-term, highly informative clinical evidence -> genomic gain contracts toward the decision boundary.

The manuscript should test this rather than merely assert it.

## Theoretical predictions to test directly

The theorem should generate empirical tests, not only post-hoc interpretation.

For clinical likelihood ratio $L_C$, conditional genomic likelihood ratio $L_G$, genomic bounds $m \leq L_G \leq M$, and decision threshold $\lambda$, the theorem predicts that a genomic contribution can alter the optimal decision only when

$$
\frac{\lambda}{M} \leq L_C \leq \frac{\lambda}{m}.
$$

Planned tests:

- [ ] Estimate or approximate the relevant genomic likelihood-ratio range for MUC5B/IPF and APOE/ADRD.
- [ ] Map ZeBRA probabilities to calibrated odds / likelihood-ratio scale.
- [ ] Determine whether observed genomic rescues are concentrated in the theoretically permitted band.
- [ ] Measure the fraction of patients lying in that band as a function of operating threshold.
- [ ] Show that incremental genomic utility contracts as the clinical score moves farther from the decision boundary.
- [ ] At increasingly stringent specificity, test whether the relative genomic contribution approaches zero compared with the clinical likelihood ratio.

## Predictor-strength analysis

The theory is not simply "EHR beats genomics." It predicts a dependence on **clinical predictor strength**.

Where possible, deliberately compare clinical models of increasing quality:

- [ ] simpler diagnosis-count / PheRS-like baseline;
- [ ] standard multivariable EHR model;
- [ ] ZeBRA.

Prediction: genomic augmentation should appear more valuable for weaker clinical predictors and progressively less valuable as the clinical predictor becomes more discriminative and its likelihood-ratio margins become larger.

This experiment would directly demonstrate why the conclusion differs from studies using more modest EHR risk scores.

## Figures that would make the paper

- [ ] ZeBRA vs genomic vs combined AUC distributions for both diseases.
- [ ] Local genomic odds ratio / information as a function of ZeBRA percentile.
- [ ] Sensitivity improvement at matched FPR across operating thresholds.
- [ ] LR+ and LR- curves for ZeBRA and genomic-boundary hybrid.
- [ ] Fraction of patients for whom genomics changes the decision vs clinical threshold.
- [ ] Theoretical boundary band overlaid on empirical genomic-rescue regions.
- [ ] Prediction-horizon curves showing genomic incremental value vs horizon.
- [ ] Cross-disease comparison: MUC5B/IPF and APOE/ADRD.
- [ ] External replication panels for Colorado vs All of Us.

## Manuscript-level positioning

Potential concise framing for the Introduction / Discussion:

> Genomic risk scores quantify inherited susceptibility, whereas longitudinal clinical histories record the evolving realization of disease. These information channels need not remain equally useful throughout the course of prediction. We show that for non-Mendelian disease, once a near-term clinical predictor enters a sufficiently informative likelihood-ratio regime, bounded genomic evidence provides little global incremental predictive value and can materially affect decisions only near the clinical operating boundary.

Potential headline sentence:

> **The genome predicts susceptibility; a sufficiently informative clinical history predicts disease. Once the latter becomes decisive, the former adds little.**

Alternative concise statement:

> **Genomic complementarity is a regime, not a universal property: it contracts as longitudinal clinical prediction becomes sufficiently accurate.**

## Evidence standard before high-profile submission

Minimum strong paper:

- [ ] theorem and complete proof;
- [x] IPF/ILD Colorado analysis;
- [ ] ADRD Colorado replication;
- [ ] strong genomic comparator for each disease (MUC5B + IPF PRS; APOE + ADRD PRS);
- [ ] fully repeated / cross-fitted global and boundary analyses;
- [ ] prespecified boundary analysis avoiding post-hoc region selection;
- [ ] explicit comparison to Detrois et al.;
- [ ] horizon-dependent analysis if technically feasible.

Flagship version:

- [ ] all of the above;
- [ ] independent All of Us replication for IPF/ILD;
- [ ] independent All of Us replication for ADRD;
- [ ] same boundary-localization geometry in both diseases and both cohorts;
- [ ] quantitative agreement between theorem-predicted boundary bands and observed genomic utility;
- [ ] demonstration that genomic incremental value diminishes as clinical predictor quality / near-term informativeness increases.

## Reference

Detrois, K. E., Hartonen, T., Teder-Laving, M., et al. Cross-biobank generalizability and accuracy of electronic health record-based predictors compared to polygenic scores. *Nature Genetics* **57**, 2136–2145 (2025). https://doi.org/10.1038/s41588-025-02298-9
