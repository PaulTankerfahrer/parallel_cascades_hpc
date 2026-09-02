#!/usr/bin/env python3
"""
analyze_fixed_energy.py -- Uebersicht ueber das Fixed-Energy-Ensemble.

Fittet fuer jede Energie und jede Clusterdefinition (Schwelle x Linkradius)
ein Potenzgesetz nach Clauset und stellt die Ergebnisse nebeneinander. Ohne
Bootstrap/GoF -- das ist die schnelle Uebersicht, mit der man die
Clusterdefinition auswaehlt. Die vertiefte Analyse macht danach
clauset_refit.py auf der ausgewaehlten Kombination.

Aufruf:
  .venv-powerlaw/bin/python scripts/python/statistik/analyze_fixed_energy.py \
      results/statistik/ensemble_fest --outdir results/statistik/fit_fest
"""
import argparse, glob, os, re, sys
import numpy as np
import powerlaw
from scipy.optimize import minimize_scalar
from scipy.special import zeta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clauset_refit import load, fit_of, pl_stats  # noqa: E402


def mle_unbounded(d, xmin):
    """Diskrete MLE ohne Parametergrenze.

    Das powerlaw-Paket sucht alpha nur in [0, 3] und liefert bei steileren
    Verteilungen stumpf 3.000 zurueck. Dieser Kreuzcheck deckt das auf.
    """
    d = d[d >= xmin]
    n = len(d)
    s = np.log(d).sum()
    r = minimize_scalar(lambda a: n * np.log(zeta(a, xmin)) + a * s,
                        bounds=(1.001, 30.0), method="bounded",
                        options=dict(xatol=1e-6))
    return float(r.x)


def combo_label(tag):
    m = re.match(r"clusters_thr(\d+)_link(\d+)", tag)
    thr, link = int(m.group(1)) / 100, int(m.group(2)) / 100
    return f"disp>{thr:.1f} / link {link:.1f}", thr, link


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ensdir")
    ap.add_argument("--outdir", default="results/statistik/fit_fest")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    edirs = sorted(glob.glob(os.path.join(a.ensdir, "E*")))
    rows = []
    for ed in edirs:
        E = int(os.path.basename(ed)[1:])
        meta = np.genfromtxt(os.path.join(ed, "run_meta.csv"), delimiter=",",
                             skip_header=1)
        nmoved = meta[:, 5]
        for fn in sorted(glob.glob(os.path.join(ed, "clusters_*.csv"))):
            tag = os.path.basename(fn)[:-4]
            lbl, thr, link = combo_label(tag)
            d = load(fn)
            if len(d) < 30:
                rows.append(dict(E=E, combo=lbl, thr=thr, link=link, n=len(d),
                                 per_run=len(d) / len(meta), skip=True))
                continue
            fit = fit_of(d)
            alpha, xmin, sigma, D, ntail = pl_stats(fit)
            alpha_free = mle_unbounded(d, xmin)
            clamped = abs(alpha - 3.0) < 1e-3 and alpha_free > 3.0 + 1e-3
            try:
                R_ln, p_ln = fit.distribution_compare(
                    "power_law", "lognormal", nested=False, normalized_ratio=True)
            except Exception:
                R_ln, p_ln = float("nan"), float("nan")
            try:
                R_tp, p_tp = fit.distribution_compare(
                    "power_law", "truncated_power_law", nested=True)
            except Exception:
                R_tp, p_tp = float("nan"), float("nan")
            rows.append(dict(E=E, combo=lbl, thr=thr, link=link, n=len(d),
                             per_run=len(d) / len(meta), nmoved=nmoved.mean(),
                             alpha=alpha, sigma=sigma, xmin=xmin, ntail=ntail,
                             D=D, xmax=int(d.max()),
                             R_ln=R_ln, p_ln=p_ln, R_tp=R_tp, p_tp=p_tp,
                             alpha_free=alpha_free, clamped=clamped,
                             tp_bad=(R_tp > 0),
                             skip=False))

    md = os.path.join(a.outdir, "fixed_energy_uebersicht.md")
    with open(md, "w") as f:
        f.write("# Fixed-Energy-Ensemble: Uebersicht\n\n")
        f.write("Feste PKA-Energie, variiert werden nur Einschlagort, Winkel "
                "und Seed. Die Einschlagparameter haengen nur von der Task-ID "
                "ab, sind also bei allen drei Energien identisch (gepaarter "
                "Vergleich).\n\n")
        f.write("`R` < 0 heisst: die Alternative beschreibt die Daten besser "
                "als das Potenzgesetz.\n\n")
        f.write("\u26a0 = das powerlaw-Paket sucht `alpha` nur in [0, 3] und "
                "haette 3.000 gemeldet; angegeben ist die freie MLE.\n"
                "`n. k.` = der Cutoff-Fit ist nicht konvergiert (das genestete "
                "Modell kann nie schlechter sein als das Potenzgesetz, ein "
                "positives R ist daher ein numerischer Fehlschlag).\n\n")
        for E in sorted({r["E"] for r in rows}):
            sub = [r for r in rows if r["E"] == E]
            nm = next((r.get("nmoved") for r in sub if r.get("nmoved")), float("nan"))
            f.write(f"\n## E = {E}  (im Mittel {nm:.0f} verlagerte Atome pro Lauf)\n\n")
            f.write("| Clusterdefinition | Cluster | pro Lauf | groesster | S | sigma | "
                    "x_min | n_tail | KS-D | R Lognormal | R Cutoff |\n")
            f.write("|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n")
            _ = None
            for r in sub:
                if r["skip"]:
                    f.write(f"| {r['combo']} | {r['n']} | {r['per_run']:.1f} | "
                            "- | - | - | - | - | - | - | - |\n")
                    continue
                a = (f"{r['alpha_free']:.3f} \u26a0" if r["clamped"]
                     else f"{r['alpha']:.3f}")
                rtp = "n.\u00a0k." if r["tp_bad"] else f"{r['R_tp']:+.1f}"
                f.write(f"| {r['combo']} | {r['n']} | {r['per_run']:.1f} | "
                        f"{r['xmax']} | {a} | {r['sigma']:.3f} | "
                        f"{r['xmin']:.0f} | {r['ntail']} | {r['D']:.4f} | "
                        f"{r['R_ln']:+.1f} | {rtp} |\n")
        f.write("\n\n## Dieselbe Clusterdefinition ueber die Energien\n\n")
        combos = []
        for r in rows:
            if r["combo"] not in combos:
                combos.append(r["combo"])
        for c in combos:
            sub = [r for r in rows if r["combo"] == c and not r["skip"]]
            if len(sub) < 2:
                continue
            al = [r["alpha_free"] if r["clamped"] else r["alpha"] for r in sub]
            f.write(f"- **{c}**: S = " +
                    ", ".join(f"{(r['alpha_free'] if r['clamped'] else r['alpha']):.2f}"
                              f" (E={r['E']})" for r in sub) +
                    f"  -- Spannweite {max(al)-min(al):.2f}\n")
    print("geschrieben:", md)
    for r in rows:
        if not r["skip"]:
            mark = " (geklemmt, frei: %.3f)" % r["alpha_free"] if r["clamped"] else ""
            print(f"E={r['E']:5d}  {r['combo']:22s}  n={r['n']:5d}  "
                  f"S={r['alpha']:.3f}  x_min={r['xmin']:.0f}  D={r['D']:.4f}{mark}")


if __name__ == "__main__":
    main()
