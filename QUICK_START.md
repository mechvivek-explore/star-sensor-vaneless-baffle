# Quick Start — New to this work?

## What problem is being studied?
A star sensor must reject strong off-axis light while preserving its required field of view. This study asks whether a continuous vane-free reflective wall can control the path of unwanted reflected rays.

## What was compared?
Five axisymmetric baffle profiles (G1–G5): straight conical, outward-curved quadratic, and inward-curved quadratic profiles.

## What is the numerical method?
A 3-D non-sequential ray-tracing (NSRT) model follows each incident ray through repeated ideal-specular wall reflections until it reaches the pupil, escapes through the entrance, or reaches the reflection cap.

## What should I look at first?
1. `figures/01_baffle_geometries.png` — understand G1–G5.
2. `figures/02_pupil_coupling_vs_angle.png` — see the principal angular result.
3. `figures/03_representative_ray_paths.png` — understand the redirection mechanism.
4. `data/key_results.csv` — see the compact numerical comparison.
5. `colab/NSRT_Quick_Start_Colab.ipynb` — explore the data in a browser.

## Main observation
For the investigated geometries, outward-curved G3/G4 reduce geometrical pupil coupling relative to same-entrance-area straight references, while inward-curved G5 shows the opposite trend. Thus the result depends on wall-profile shape, not curvature alone.

These are ideal-specular geometrical results, not PST/NPST or detector-radiometric predictions.
