# Data guide

The CSV files in this directory are the compact, publication-oriented numerical
resources for the study.

- `geometry_definitions.csv`: frozen G1–G5 geometry parameters.
- `full_nsrt_results.csv`: complete 0–85° grouped NSRT result set (5 geometries × 122 sampled angles).
- `key_results.csv`: compact values used in the manuscript comparisons.

The numerical results are from ideal-specular geometrical NSRT. They should not be
interpreted as coating BRDF, PST/NPST, detector irradiance, or radiometric performance.

For `full_nsrt_results.csv`, ray counts satisfy:
`P + E + T50plus + X = launched`
for every archived case, with `X = 0`.
