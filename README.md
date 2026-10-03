# Continuous Vaneless Reflective Baffle for a Star Sensor

**3-D non-sequential ray tracing (NSRT) · reusable code · numerical data · guided Colab exploration**

Frozen reproducibility and research-use archive accompanying:

**“तारा संवेदक हेतु सतत वैन-रहित परावर्तक बैफल: एक वैकल्पिक अभिकल्पना”**  
**Vivek Shukla — U R Rao Satellite Centre (URSC), ISRO, Bengaluru, India**

## ⚡ Need the code immediately?

**Core reusable NSRT engine → [code/nsrt_3d.py](code/nsrt_3d.py)**  
Production driver → [code/run_analysis.py](code/run_analysis.py)  
Verification → [code/verify_results.py](code/verify_results.py)  
Model assumptions/settings → [METHOD_AND_SETTINGS.md](METHOD_AND_SETTINGS.md)

**No local Python setup?** Open [colab/NSRT_Quick_Start_Colab.ipynb](colab/NSRT_Quick_Start_Colab.ipynb) in Google Colab and explore the data/model from a browser.

## 🔎 New to this research area?

Start with **[QUICK_START.md](QUICK_START.md)**. It guides you through the physical problem, G1–G5 geometries, main result, ray-redirection mechanism, data, and finally the NSRT code.

## Main observation

For the investigated geometries, outward-curved G3/G4 profiles reduce geometrical pupil coupling relative to same-entrance-area straight references, while inward-curved G5 shows the opposite trend. The result therefore depends on **wall-profile shape**, not curvature alone.

## Repository map

- **code/** — reusable 3-D NSRT engine, production driver, verification.
- **data/** — geometry definitions, complete 0–85° publication-analysis dataset, key comparisons.
- **figures/** — three visual entry points to the study.
- **colab/** — browser-based guided exploration.
- **METHOD_AND_SETTINGS.md** — assumptions, run settings, ray-fate definitions.
- **QUICK_START.md** — guided route for newcomers.

## Data for further research

- `data/geometry_definitions.csv` — G1–G5 parameters.
- `data/full_nsrt_results.csv` — consolidated 0–85° NSRT results.
- `data/key_results.csv` — compact comparisons used in the paper.

The calculations assume ideal specular reflection and report geometrical pupil coupling. They are not PST/NPST, coating-BRDF, detector-irradiance, or radiometric predictions.

## Verification

Run:

`python code/verify_results.py`

Expected result:

`PASS: conservation and key manuscript values.`

## Frozen release and reuse note

This repository is intended to be a **frozen publication archive**. The published version is not intended for routine editing. Researchers may download or clone the archive for independent review, verification and further development at their end.

No reuse licence is asserted here. Any licence or public-reuse permission should be added only after confirming the applicable organisational/publication/IP requirements.
