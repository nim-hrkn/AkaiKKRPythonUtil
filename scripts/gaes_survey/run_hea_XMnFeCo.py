"""GAES (Method 2) for X-Mn-Fe-Co fcc HEAs, X from the shallow-core-level table (no noble gases, no reconf stoppers)."""
import json, os, sys, time, traceback
from pyakaikkr.gaes import Gaes, Layout, SiteComposition, make_single_site_param
EXE = "/home/kino/kino/Claude/Project/AKAIKKR/AkaiKKRprogram.2022.0721.ifort/akaikkr/specx"
X = ["Ga", "As", "Se", "Rb", "In", "Sb", "Te", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Yb", "Lu", "Tl", "Bi"]
os.environ["OMP_NUM_THREADS"] = "20"
summary = {}
for x in X:
    name = x + "MnFeCo"
    comp = SiteComposition.from_elements([x, "Mn", "Fe", "Co"])
    prefix = os.path.join("RUN", name + "_fcc")
    t0 = time.time()
    try:
        params = {"fcc": make_single_site_param(comp, "fcc", type_name="HEA")}
        g = Gaes(EXE, Layout(prefix, version=2), ewidth_init=1.2, ewidth_dos=3.0, ref=0.75, method=2,
                 dosth=2e-2, dosth2=1e-3, min_ewidth=1.0, max_ewidth=2.0, max_pm_iter=2, with_j=False)
        res = g.run(name, params)
        summary[name] = dict(status=res.status, ewidth=res.ewidth_final, tried=res.ewidth_tried, gap=res.gap_used,
                             flags=[j.flag for j in res.judgements], converged=res.converged, message=res.message,
                             seconds=round(time.time() - t0))
    except Exception as e:
        summary[name] = dict(status="exception", message=repr(e), seconds=round(time.time() - t0))
        traceback.print_exc()
    print(name, summary[name], flush=True)
    json.dump(summary, open("summary.json", "w"), indent=1)
