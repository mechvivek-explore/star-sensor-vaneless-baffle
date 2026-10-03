"""Check conservation and key manuscript values in the public dataset."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
df=pd.read_csv(ROOT/"data/full_nsrt_results.csv")
assert set(df.geometry_id)=={"G1","G2","G3","G4","G5"}
assert df.groupby("geometry_id").angle_deg.nunique().eq(122).all()
assert (df["P"]+df["E"]+df["T50plus"]+df["X"]==df["launched"]).all()
assert (df["X"]==0).all()
checks={
("G1",24):2693,("G3",24):1748,("G2",24):2678,("G4",24):1449,("G5",24):2929,
("G1",28):709,("G3",28):254,("G2",28):674,("G4",28):107,("G5",28):884,
}
for (g,a),expected in checks.items():
    got=int(df.loc[(df.geometry_id==g)&(df.angle_deg==a),"P"].iloc[0])
    assert got==expected,(g,a,got,expected)
print("PASS: conservation and key manuscript values.")
