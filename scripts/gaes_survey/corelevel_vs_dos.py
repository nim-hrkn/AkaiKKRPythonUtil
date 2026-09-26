"""core levels printed in out_go.log vs. semicore peaks of the total DOS, over the 2019 RUN directories.

Usage (in fukushima_HEA_run_exprlattice/production_run/run0): python corelevel_vs_dos.py
Writes corelevel_vs_dos.json to GAES_SCRATCH (default .) and prints the per-(element, level) statistics.
See docs/core_levels_vs_dos.md.
"""
import os, re, sys, json, collections
import numpy as np
from multiprocessing import Pool
S = os.environ.get("GAES_SCRATCH", ".")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from scan_semicore import read_total_dos
RUN="RUN"
def analyse(d):
    try:
        txt=open(os.path.join(RUN,d,"out_go.log")).read()
        if "sbtime report" not in txt: return None
        ef=float(re.findall(r"\n\s+ef=\s*([-\d.]+)",txt)[-1])
        r=read_total_dos(os.path.join(RUN,d,"out_dos.log"))
        if r is None: return None
        e,s=r
        # local maxima of the DOS below the valence band bottom region: take all local maxima with height>1
        pk=[(e[i],s[i]) for i in range(1,len(e)-1) if s[i]>s[i-1] and s[i]>=s[i+1] and s[i]>1.0 and e[i]<-0.85]
        out=[]
        for b in re.split(r"\*\*\* type-",txt)[1:]:
            m=re.match(r"\S+\s+(\S+)\s+\(z=\s*([\d.]+)\)",b)
            if not m: continue
            el=m.group(1)
            for v,l,st in re.findall(r"(-?\d+\.\d+) Ry\((\d[spdf])\)(\*?)",b):
                er=float(v)-ef
                if er<-2.3 or er>0: continue
                near=min(pk,key=lambda p:abs(p[0]-er)) if pk else None
                out.append((el,l,bool(st),round(er,3),round(near[0],3) if near else None,round(near[1],1) if near else None))
        return d,out
    except Exception:
        return None
if __name__=="__main__":
    dirs=sorted(x for x in os.listdir(RUN) if x.endswith("polytyp_fcc,pm_000"))
    with Pool(20) as p: res=[r for r in p.map(analyse,dirs,chunksize=40) if r]
    json.dump(res, open(os.path.join(S, "corelevel_vs_dos.json"), "w"))
    # statistics per (element, level, valence flag)
    stat=collections.defaultdict(list)
    for d,out in res:
        for el,l,st,er,pe,ph in out:
            if pe is not None: stat[(el,l,st)].append((er,pe-er,ph))
    print("dirs analysed:",len(res))
    print("%-4s %-3s %-8s %6s %10s %10s %10s %8s"%("el","lvl","treated","n","level_med","peak-lvl","|d|>0.1","peak_h"))
    for key in sorted(stat,key=lambda k:(k[0],k[1])):
        a=np.array(stat[key]); 
        if len(a)<20: continue
        print("%-4s %-3s %-8s %6d %10.2f %10.3f %9.0f%% %8.0f"%(key[0],key[1],"valence*" if key[2] else "core",len(a),np.median(a[:,0]),np.median(a[:,1]),100*np.mean(np.abs(a[:,1])>0.1),np.median(a[:,2])))
