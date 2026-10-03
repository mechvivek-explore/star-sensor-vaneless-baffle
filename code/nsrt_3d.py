"""3-D non-sequential ideal-specular ray tracer for axisymmetric star-sensor baffles.

z=0: pupil plane; z=L: entrance plane. Wall is an ordered (z,r) polyline revolved
about z. Every propagation step evaluates ALL wall/aperture candidates and selects
the nearest positive physical event. Ray fractions are geometrical, not PST/NPST.
"""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional
import math, csv
import numpy as np

EPS=1e-9; TOL=1e-8; EVENT_TIE_TOL=1e-8; MAX_REFLECTIONS=50

class Fate(Enum): PUPIL=auto(); ESCAPE=auto(); TRAPPED=auto(); ERROR=auto()
@dataclass(frozen=True)
class Hit:
    t: float; kind: str; point: np.ndarray; normal: Optional[np.ndarray]=None; segment: Optional[int]=None
@dataclass(frozen=True)
class RayResult:
    fate:Fate; reflections:int; point:np.ndarray; path:tuple; segment_history:tuple
@dataclass
class AngleResult:
    angle_deg:float; launched:int; P:np.ndarray; E:np.ndarray; T:int; X:int
    @property
    def pupil_total(self): return int(self.P.sum())
    @property
    def escape_total(self): return int(self.E.sum())
    @property
    def pupil_fraction(self): return self.pupil_total/self.launched
    @property
    def escape_fraction(self): return self.escape_total/self.launched
    @property
    def checksum(self): return self.pupil_total+self.escape_total+self.T+self.X

def unit(v):
    v=np.asarray(v,float); n=float(np.linalg.norm(v))
    if not np.isfinite(n) or n<=0: raise ValueError('zero/non-finite vector')
    return v/n

def reflect(s,n):
    s=unit(s); n=unit(n); return unit(s-2*np.dot(s,n)*n)

def sunflower_disk(n,radius,phase=0.0):
    if n<=0 or radius<=0: raise ValueError
    i=np.arange(n,dtype=float); rho=radius*np.sqrt((i+0.5)/n)
    phi=phase+i*math.pi*(3-math.sqrt(5))
    return np.column_stack((rho*np.cos(phi),rho*np.sin(phi)))

def source_direction(theta_deg,azimuth_deg=0.0):
    th=math.radians(theta_deg); ph=math.radians(azimuth_deg)
    return unit([math.sin(th)*math.cos(ph),math.sin(th)*math.sin(ph),-math.cos(th)])

class AxisymmetricBaffle:
    def __init__(self,profile_zr,pupil_radius,entrance_z,entrance_radius):
        p=np.asarray(list(profile_zr),float)
        if p.ndim!=2 or p.shape[1]!=2 or len(p)<2 or not np.all(np.isfinite(p)): raise ValueError('bad profile')
        if np.any(p[:,1]<=0): raise ValueError('r must be >0')
        self.profile=p; self.Rp=float(pupil_radius); self.L=float(entrance_z); self.Re=float(entrance_radius)
    @staticmethod
    def cone(Rp,L,Re,n_segments=1):
        z=np.linspace(0,L,n_segments+1); r=Rp+(Re-Rp)*z/L
        return AxisymmetricBaffle(zip(z,r),Rp,L,Re)

    def wall_candidates(self,o,d):
        hits=[]; o=np.asarray(o,float); d=np.asarray(d,float)
        ox,oy,oz=map(float,o); dx,dy,dz=map(float,d)
        z0=self.profile[:-1,0]; r0=self.profile[:-1,1]; z1=self.profile[1:,0]; r1=self.profile[1:,1]; ds=z1-z0
        const=np.abs(ds)<=TOL
        if np.any(const):
            for j in np.flatnonzero(const):
                if abs(dz)<=TOL: continue
                t=(z0[j]-oz)/dz
                if t<=EPS: continue
                q=o+t*d; rq=math.hypot(q[0],q[1]); lo,hi=sorted((r0[j],r1[j]))
                if lo-TOL<=rq<=hi+TOL: hits.append(Hit(float(t),'wall',q,np.array([0.,0.,1.]),int(j)))
        idx=np.flatnonzero(~const)
        if idx.size:
            zz0=z0[idx]; zz1=z1[idx]; rr0=r0[idx]; rr1=r1[idx]; dss=ds[idx]
            a=(rr1-rr0)/dss; bb=rr0-a*zz0; c0=a*oz+bb; c1=a*dz
            A=dx*dx+dy*dy-c1*c1; B=2*(ox*dx+oy*dy-c0*c1); C=ox*ox+oy*oy-c0*c0
            linear=np.abs(A)<=1e-14; candidates=[]
            validlin=linear & (np.abs(B)>1e-14)
            if np.any(validlin):
                for k,t in zip(np.flatnonzero(validlin),(-C[validlin]/B[validlin])): candidates.append((k,float(t)))
            disc=B*B-4*A*C; validq=(~linear) & (disc>=-1e-12)
            if np.any(validq):
                ks=np.flatnonzero(validq); sd=np.sqrt(np.maximum(0.,disc[validq]))
                for k,t1,t2 in zip(ks,(-B[validq]-sd)/(2*A[validq]),(-B[validq]+sd)/(2*A[validq])):
                    candidates.append((k,float(t1))); candidates.append((k,float(t2)))
            for k,t in candidates:
                if not np.isfinite(t) or t<=EPS: continue
                q=o+t*d; lo=min(zz0[k],zz1[k]); hi=max(zz0[k],zz1[k])
                if not(lo-TOL<=q[2]<=hi+TOL): continue
                rq=math.hypot(q[0],q[1]); er=a[k]*q[2]+bb[k]
                if er<=0 or abs(rq-er)>5e-7 or rq<=TOL: continue
                n=unit([q[0],q[1],-a[k]*er]); hits.append(Hit(float(t),'wall',q,n,int(idx[k])))
        return hits

    def aperture_candidates(self,o,d):
        hits=[]; dz=float(d[2])
        if abs(dz)<=TOL: return hits
        tp=-o[2]/dz
        if tp>EPS:
            q=o+tp*d
            if math.hypot(q[0],q[1])<=self.Rp+TOL: hits.append(Hit(float(tp),'pupil',q))
        if dz>0:
            te=(self.L-o[2])/dz
            if te>EPS:
                q=o+te*d
                if math.hypot(q[0],q[1])<=self.Re+TOL: hits.append(Hit(float(te),'entrance',q))
        return hits

    def next_event(self,o,d):
        c=self.wall_candidates(o,d)+self.aperture_candidates(o,d)
        if not c: return None
        tmin=min(h.t for h in c)
        tied=[h for h in c if h.t <= tmin + EVENT_TIE_TOL]
        priority={'pupil':0,'entrance':1,'wall':2}
        return min(tied,key=lambda h:(priority[h.kind],h.t,h.segment if h.segment is not None else -1))

    def trace(self,origin,direction,max_reflections=MAX_REFLECTIONS,keep_path=False):
        o=np.asarray(origin,float).copy(); d=unit(direction); nr=0; path=[o.copy()] if keep_path else []; hist=[]
        while True:
            h=self.next_event(o,d)
            if h is None: return RayResult(Fate.ERROR,nr,o.copy(),tuple(path),tuple(hist))
            if keep_path: path.append(h.point.copy())
            if h.kind=='pupil': return RayResult(Fate.PUPIL,nr,h.point.copy(),tuple(path),tuple(hist))
            if h.kind=='entrance': return RayResult(Fate.ESCAPE,nr,h.point.copy(),tuple(path),tuple(hist))
            if nr>=max_reflections: return RayResult(Fate.TRAPPED,nr,h.point.copy(),tuple(path),tuple(hist))
            d=reflect(d,h.normal); nr+=1; hist.append(h.segment); o=h.point+1e-8*d

    def launch_bundle(self,theta_deg,n_rays,azimuth_deg=0.0,phase=0.0,max_reflections=MAX_REFLECTIONS):
        xy=sunflower_disk(n_rays,self.Re,phase); d=source_direction(theta_deg,azimuth_deg)
        crossings=np.column_stack((xy[:,0],xy[:,1],np.full(n_rays,self.L)))
        origins=crossings+1e-7*d
        P=np.zeros(max_reflections+1,dtype=np.int64); E=P.copy(); T=X=0
        for o in origins:
            rr=self.trace(o,d,max_reflections)
            if rr.fate is Fate.PUPIL: P[rr.reflections]+=1
            elif rr.fate is Fate.ESCAPE: E[rr.reflections]+=1
            elif rr.fate is Fate.TRAPPED: T+=1
            else: X+=1
        out=AngleResult(theta_deg,n_rays,P,E,T,X)
        assert out.checksum==n_rays; assert E[0]==0
        return out

def sweep(baffle,angles_deg,n_rays,azimuth_deg=0.0,phase=0.0,max_reflections=MAX_REFLECTIONS):
    return [baffle.launch_bundle(a,n_rays,azimuth_deg,phase,max_reflections) for a in angles_deg]

def write_csv(results,filename):
    m=len(results[0].P)-1
    fields=['angle_deg','launched','pupil_total','escape_total','T','X','pupil_fraction','escape_fraction','P_C_0_to_8','P_B_9_to_50','E_A_1_to_8','E_B_9_to_50']+[f'P{i}' for i in range(m+1)]+[f'E{i}' for i in range(m+1)]
    with open(filename,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in results:
            row=dict(angle_deg=r.angle_deg,launched=r.launched,pupil_total=r.pupil_total,escape_total=r.escape_total,T=r.T,X=r.X,pupil_fraction=r.pupil_fraction,escape_fraction=r.escape_fraction,P_C_0_to_8=int(r.P[:9].sum()),P_B_9_to_50=int(r.P[9:].sum()),E_A_1_to_8=int(r.E[1:9].sum()),E_B_9_to_50=int(r.E[9:].sum()))
            row.update({f'P{i}':int(r.P[i]) for i in range(m+1)}); row.update({f'E{i}':int(r.E[i]) for i in range(m+1)}); w.writerow(row)

def run_self_tests():
    rng=np.random.default_rng(1234); ne=de=0.
    for _ in range(10000):
        s=unit(rng.normal(size=3)); n=unit(rng.normal(size=3)); r=reflect(s,n)
        ne=max(ne,abs(np.linalg.norm(r)-1)); de=max(de,abs(np.dot(r,n)+np.dot(s,n)))
    assert ne<2e-14 and de<2e-14
    b=AxisymmetricBaffle.cone(10,100,30); h=b.next_event(np.array([0.,0.,50.]),np.array([1.,0.,0.])); assert h.kind=='wall' and np.linalg.norm(h.point-[20,0,50])<1e-9
    rr=AxisymmetricBaffle.cone(10,100,40).trace(np.array([0.,0.,100-1e-7]),np.array([0.,0.,-1.])); assert rr.fate is Fate.PUPIL and rr.reflections==0
    Re=10+100*math.tan(math.radians(16)); c=AxisymmetricBaffle.cone(10,100,Re); a=c.launch_bundle(40,257); assert a.checksum==257 and a.E[0]==0
    print('SELF-TESTS PASS'); print('reflection unit error',ne); print('reflection normal-component error',de); print('analytic cone hit',h.point); print('40-deg checksum',a.checksum,'P/E/T/X',a.pupil_total,a.escape_total,a.T,a.X)

if __name__=='__main__':
    run_self_tests()
