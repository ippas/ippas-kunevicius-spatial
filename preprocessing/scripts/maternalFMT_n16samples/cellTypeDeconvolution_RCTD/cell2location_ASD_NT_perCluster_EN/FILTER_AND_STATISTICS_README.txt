CELL-TYPE SELECTION AND STATISTICAL TESTING
===========================================

There is NO universal, officially prescribed cell2location cutoff for
filtering rare cell populations before testing group differences.

To avoid arbitrary thresholds such as 1% or 5%, as well as any selection
based on p-values, this pipeline uses the following transparent approach:

1. TESTING: Include all 59 reference populations.

2. WHOLE-SECTION COMPARISON: For each population, compare mean estimated
   section-level proportions between ASD and Neurotypical. Adjust p-values
   across ALL testable populations using Benjamini-Hochberg FDR.

3. PER-CLUSTER COMPARISON: Compare each population within each cluster,
   using section-level means. The PRIMARY FDR adjustment spans ALL testable
   cluster-by-population hypotheses globally. FDR adjusted separately within
   each cluster is also saved as a secondary result.

4. FIGURES: Show the ten most abundant populations across all 16 sections,
   ranked without knowing ASD/Neurotypical group labels or p-values.
   This selection is for visualization only. Do not recalculate FDR for the
   top ten or present them as if they were the only tested populations.

5. SAMPLE PRESENCE: A cluster comparison requires at least three sections
   from EACH group containing that cluster; otherwise the test is not
   performed and the p-value remains missing (NA).

6. BIOLOGICAL REPLICATION: Multiple sections from the same mouse are not
   independent biological replicates. Section-level p-values assume
   independence and are exploratory. Confirmatory testing should account
   for animal IDs by aggregating sections by animal or fitting an appropriate
   hierarchical model.

REFERENCES

Official cell2location tutorial (including gene filtering and posterior q05;
no universal rare-cell-type test filter):
  https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html

Official mouse brain example (published reference signatures):
  https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo.html

Propeller (testing cell-type proportions across biological replicates):
  https://doi.org/10.1093/bioinformatics/btac582

NOTE: This is a documented analysis choice, not a claimed official
"standard filter" copied from the literature.
