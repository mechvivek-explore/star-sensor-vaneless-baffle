"""Reproduce the public 0–85° geometrical NSRT dataset for G1–G5.

The calculation uses the frozen RS-01/RS-01A settings described in
METHOD_AND_SETTINGS.md. Output is written inside the repository by default.

This is geometrical ideal-specular ray tracing; it is not PST/NPST or a
detector-radiometric model.
"""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path
import numpy as np
from nsrt_3d import AxisymmetricBaffle, MAX_REFLECTIONS

RP=10.0
L=100.0
N0=20000
R0=37.0
PHASE=0.37
AZIMUTH_DEG=0.0
QUAD_INTERVALS=500

PRODUCTION_ANGLES=np.array(
    [15.,24.,28.,29., *np.round(np.arange(30.0,40.0+0.0001,0.1),10), 60.,85.]
)
FULL_ANGLES=np.unique(np.concatenate((np.arange(0.,15.,1.),PRODUCTION_ANGLES)))
assert len(PRODUCTION_ANGLES)==107
assert len(FULL_ANGLES)==122

GEOMETRIES={
 'G1': dict(name='16° Straight Conical Baffle', a=0.28674538575880792, b=0.0, Re=38.674538575880788),
 'G2': dict(name='15° Straight Conical Baffle', a=0.2679491924311227, b=0.0, Re=36.794919243112268),
 'G3': dict(name='Outward-Curved Quadratic Baffle – Large Entrance', a=0.39, b=-0.0010325461424119211, Re=38.674538575880788),
 'G4': dict(name='Outward-Curved Quadratic Baffle – Minimum Entrance', a=0.39, b=-0.0012205080756887731, Re=36.794919243112268),
 'G5': dict(name='Inward-Curved Quadratic Baffle – Large Entrance', a=0.2679491924311227, b=0.00018796193327685201, Re=38.674538575880788),
}

def launch_count(Re):
    return int(round(N0*(Re/R0)**2))

def build_baffle(g):
    z=np.array([0.0,L]) if abs(g['b'])<1e-30 else np.linspace(0.0,L,QUAD_INTERVALS+1)
    r=RP+g['a']*z+g['b']*z*z
    assert abs(r[0]-RP)<1e-12 and abs(r[-1]-g['Re'])<1e-10
    return AxisymmetricBaffle(zip(z,r),RP,L,g['Re'])

def summarize(gid, angle, rr):
    return dict(
        geometry_id=gid, angle_deg=float(angle), launched=rr.launched,
        P=rr.pupil_total, P0=int(rr.P[0]), P1_8=int(rr.P[1:9].sum()),
        P9plus=int(rr.P[9:].sum()), E=rr.escape_total, E0=int(rr.E[0]),
        E1_8=int(rr.E[1:9].sum()), E9plus=int(rr.E[9:].sum()),
        T50plus=rr.T, X=rr.X
    )

def calculate(angles=FULL_ANGLES, geometry_ids=None):
    geometry_ids=geometry_ids or list(GEOMETRIES)
    rows=[]
    for gid in geometry_ids:
        g=GEOMETRIES[gid]
        baffle=build_baffle(g)
        n=launch_count(g['Re'])
        for angle in angles:
            rr=baffle.launch_bundle(float(angle),n,AZIMUTH_DEG,PHASE,MAX_REFLECTIONS)
            if rr.checksum!=n or rr.X!=0:
                raise RuntimeError(f'Failed conservation/classification: {gid}, {angle} deg')
            rows.append(summarize(gid,angle,rr))
    return rows

def write_rows(rows, output):
    fields=['geometry_id','angle_deg','launched','P','P0','P1_8','P9plus',
            'E','E0','E1_8','E9plus','T50plus','X']
    output=Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return output

def smoke_test():
    """Fast execution-path check; frozen full-density values are checked separately
    by code/verify_results.py against the archived public dataset."""
    gid='G3'
    g=GEOMETRIES[gid]
    baffle=build_baffle(g)
    n=257
    rr=baffle.launch_bundle(28.0,n,AZIMUTH_DEG,PHASE,MAX_REFLECTIONS)
    assert rr.checksum==n and rr.X==0
    assert int(rr.P.sum())+int(rr.E.sum())+rr.T+rr.X==n
    print(f"SMOKE TEST PASS: {gid}, 28 deg, {n} rays; P/E/T/X = "
          f"{rr.pupil_total}/{rr.escape_total}/{rr.T}/{rr.X}")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',default='data/reproduced_full_nsrt_results.csv')
    ap.add_argument('--smoke-test',action='store_true',
                    help='run two quick frozen-value checks instead of the full sweep')
    args=ap.parse_args()
    if args.smoke_test:
        smoke_test()
    else:
        out=write_rows(calculate(),args.output)
        print(f"Wrote {out}")

if __name__=='__main__':
    main()
