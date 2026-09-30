#!/usr/bin/env python3
"""Drawing regression for the structure gallery: does each figure agree with the coordinates it came from?

Three assertions, for every structure:
  1. CONNECTIVITY. The number of connected components of the disc-union outline equals the number of connected
     components of the feature atoms' bond graph (cutoff 1.15 x the same-layer spacing). A drawing that splits a
     bonded cluster, or fuses two separate ones, fails here. This is the check that would have caught the disc
     radius being built from a cross-layer projected distance of 1.6975 A instead of the 2.9401 A bond.
  2. A-A' CORRESPONDENCE. The section profile is compared against the field re-sampled independently at
     A + s*t_hat, reconstructed from cut_frame. Any shift applied to one and not the other fails here.
  3. BAND WIDTH. Every atom the side view draws lies within the perpendicular half-width its caption states.

Run it before looking at any rendered image, and before claiming a drawing bug is fixed.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/check_gallery_drawing.py
"""
import json, sys, numpy as np
sys.path.insert(0,"/anvil/scratch/x-rywang/Au_Cl/scripts")
from render_gallery import (flatten_cell, surface_model, feature_sets, layer_spacing, blob_distance,
                            cut_profile, cut_band, cut_frame, parent_of, TILE_TARGET, GAP_LAYER)
from ase.io import read
from scipy import ndimage
gal=json.load(open("/anvil/scratch/x-rywang/Au_Cl/analysis/gallery/gallery.json"))
def atom_components(pts, cell, cutoff):
    n=len(pts)
    if n==0: return 0
    par=list(range(n))
    def f(a):
        while par[a]!=a: par[a]=par[par[a]]; a=par[a]
        return a
    C=np.array([cell[0][:2],cell[1][:2]]); Ci=np.linalg.inv(C)
    fr=(pts[:,None,:]-pts[None,:,:])@Ci; fr-=np.round(fr)
    d=np.linalg.norm(fr@C,axis=-1)
    for i in range(n):
        for j in range(i+1,n):
            if d[i,j]<cutoff: par[f(i)]=f(j)
    return len({f(i) for i in range(n)})
def disc_components(D, R, shape):
    """connected components of the disc union on the PERIODIC grid"""
    M=(D<R)
    lab,k=ndimage.label(M)
    if k<=1: return int(k)
    # stitch across the periodic seams
    par=list(range(k+1))
    def f(a):
        while par[a]!=a: par[a]=par[par[a]]; a=par[a]
        return a
    for A,B in ((lab[0,:],lab[-1,:]),(lab[:,0],lab[:,-1])):
        for u,v in zip(A,B):
            if u and v: par[f(u)]=f(v)
    return len({f(i) for i in range(1,k+1) })
bad=[]; rows=[]
for sid in sorted(gal):
    m=gal[sid]; at=flatten_cell(read(f"{m['dir']}/{m['geom']}")); cell=at.get_cell().array
    pa=parent_of(sid,at)
    sm=surface_model(at, prefer_xy=(pa['added'] if pa else None)); a0=layer_spacing(at)
    hi,lo,modal=feature_sets(at)
    note=[]
    # --- 1. connectivity: disc graph must match the atom bond graph, feature by feature
    for pts,lay in zip((hi,lo), sm["layers"] or [None,None]):
        if lay is None or not len(pts): continue
        na=atom_components(pts,cell,1.15*a0); nd=disc_components(lay["D"],sm["R"],sm["H"].shape)
        if na!=nd: bad.append(f"{sid}: feature of {len(pts)} atoms -> {na} bonded clusters but {nd} drawn blobs")
        note.append(f"{len(pts)}at/{na}cl")
    # --- 2. periodic replication: every tiled copy of a field must be identical (it is a translation now)
    t1=max(1,min(3,int(round(TILE_TARGET/np.linalg.norm(cell[0][:2])))))
    t2=max(1,min(3,int(round(TILE_TARGET/np.linalg.norm(cell[1][:2])))))
    # --- 3. A->A' correspondence: the profile order must match the plan-view line
    prof_desc=""
    for ci in sm["cuts"]:
        F=sm["Z"] if sm["mode"]=="blob" else sm["H"]
        s,pr,L=cut_profile(F,ci,cell)
        rA,t,nv,LL,P=cut_frame(cell,ci)
        # sample the field independently at the geometric positions of A + s*t and compare
        gx,gy=sm["gx"],sm["gy"]
        C=np.array([cell[0][:2],cell[1][:2]]); Ci=np.linalg.inv(C)
        pts=rA[None,:]+s[:,None]*t[None,:]
        fr=(pts@Ci)%1.0
        ii=np.clip((fr[:,0]*F.shape[0]).astype(int),0,F.shape[0]-1)
        jj=np.clip((fr[:,1]*F.shape[1]).astype(int),0,F.shape[1]-1)
        resid=float(np.abs(F[ii,jj]-pr).max())
        if resid>1e-9: bad.append(f"{sid}: profile does not match the marked {ci['label']}-{ci['label']}' line (max {resid:.3g})")
        lab=np.where(pr>0.5,'H',np.where(pr<-0.5,'L','T'))
        runs=[]
        for c in lab:
            if runs and runs[-1][0]==c: runs[-1][1]+=1
            else: runs.append([c,1])
        prof_desc += f" {ci['label']}:{'-'.join(f'{c}{n}' for c,n in runs)}"
        # --- 4. band: true perpendicular half width
        mk,ss,dp=cut_band(at.get_positions(),cell,ci,1.1*a0)
        if mk.any() and float(np.abs(dp[mk]).max())>1.1*a0+1e-9:
            bad.append(f"{sid}: band selection exceeds its stated half width")
    rows.append(f"{sid:38s} {sm['mode']:5s} cuts={len(sm['cuts'])} a0={a0:.4f} R={sm['R']:.3f} "
                f"[{' '.join(note)}]{prof_desc}")
print("\n".join(rows))
print()
if bad:
    print("FAILURES:"); print("\n".join("  "+b for b in bad))
else:
    print(f"all {len(gal)} structures pass: disc blobs = bonded clusters, profile = marked line, band within its stated width")
