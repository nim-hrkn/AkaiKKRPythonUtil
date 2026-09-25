"""GAES (Method 2, ewidth_init 1.2, [1.0, 2.0]) on tests/structure/small_primitive_structures (akaikkr build, nmag)."""
import glob, json, os, time, traceback, warnings
warnings.filterwarnings("ignore")
from pyakaikkr.Cif2Kkr import ak_cif2kkrparam
from pyakaikkr.gaes import Gaes, Layout
EXE = "/home/kino/kino/Claude/Project/AKAIKKR/AkaiKKRprogram.2022.0721.ifort/akaikkr/specx"
SRC = "/home/kino/kino/Claude/Project/AKAIKKR/AkaiKKRPythonUtil/tests/structure/small_primitive_structures"
os.environ["OMP_NUM_THREADS"] = "20"
summary = {}
for cif in sorted(glob.glob(SRC + "/*/*_primitive.cif")):
    name = os.path.basename(os.path.dirname(cif)); t0 = time.time()
    try:
        p = ak_cif2kkrparam(cif)
        p = p[0] if isinstance(p, tuple) else p
        p.update({"go": "go", "potentialfile": "pot.dat", "magtyp": "nmag", "sdftyp": "pbe", "reltyp": "sra", "bzqlty": 8,
                  "record": "2nd", "outtyp": "update", "edelt": 1e-4, "ewidth": 1.2, "pmix": 0.005, "maxitr": 500})
        g = Gaes(EXE, Layout(os.path.join("RUN", name), version=2), ewidth_init=1.2, ewidth_dos=3.0, ref=0.75, method=2,
                 dosth=2e-2, dosth2=1e-3, min_ewidth=1.0, max_ewidth=2.0, max_pm_iter=2, with_j=False)
        res = g.run(name, {"prim": p})
        summary[name] = dict(status=res.status, ewidth=res.ewidth_final, tried=res.ewidth_tried, gap=res.gap_used,
                             flags=[j.flag for j in res.judgements], converged=res.converged, message=res.message,
                             reasons=[j.reasons for j in res.judgements], seconds=round(time.time() - t0),
                             brvtyp=p["brvtyp"], natm=p["natm"], types=p["type"])
    except Exception as e:
        summary[name] = dict(status="exception", message=repr(e), seconds=round(time.time() - t0)); traceback.print_exc()
    print(name, summary[name], flush=True)
    json.dump(summary, open("summary.json", "w"), indent=1)
