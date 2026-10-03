# Analysis settings

- Pupil radius: 10 mm
- Baffle length: 100 mm
- Full field of view: 30 deg (±15 deg)
- Profile: `r(z) = 10 + a z + b z^2`, 0 <= z <= 100 mm
- Same incident ray density used for geometry comparison
- G1/G3/G5: entrance radius 38.6745 mm; 21,851 rays
- G2/G4: entrance radius 36.7949 mm; 19,779 rays
- Ideal specular wall reflection
- Maximum wall reflections: 50
- Ray fates: pupil (`P_n`), entrance escape (`E_n`), unresolved after cap (`T50+`), numerical/unclassified (`X`)
- Production phase: 0.37 rad
- `data/full_nsrt_results.csv` contains the consolidated 0–85 deg publication-analysis dataset.

The results are geometrical pupil-coupling results, not PST/NPST or detector-radiometric predictions.

## Reproduction

`code/run_analysis.py` reproduces the public 0–85° grouped NSRT dataset using the same geometry, ray density, phase, ideal-specular reflection law and reflection cap.

For a quick installation/code check:

`python code/run_analysis.py --smoke-test`

For the complete sweep:

`python code/run_analysis.py`

The full sweep is computationally much heavier than the smoke test.
