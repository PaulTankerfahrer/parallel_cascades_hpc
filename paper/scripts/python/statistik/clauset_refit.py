#!/usr/bin/env python3
"""
clauset_refit.py -- rigoroser Potenzgesetz-Fit nach Clauset, Shalizi & Newman
(2009) mit dem Python-Paket `powerlaw` (Alstott et al. 2014).

Liefert pro Datensatz:
  1) x_min via KS-Minimierung (uneingeschraenkter Scan ueber alle Kandidaten)
  2) S (= alpha) via diskretem MLE + Standardfehler sigma = (S-1)/sqrt(n_tail)
  3) Bootstrap-Konfidenzintervall (x_min wird in jedem Resample neu geschaetzt)
  4) Goodness-of-Fit-p-Wert (semi-parametrischer Bootstrap, Clauset Abschn. 4)
  5) Likelihood-Ratio-Tests gegen Lognormal, Exponential, Potenzgesetz mit
     exponentiellem Cutoff und Stretched Exponential (Vuong-Test, Clauset Abschn. 5)
  6) Kontrollfit bei festem x_min (Vergleichbarkeit mit dem Bericht)

Aufruf:
  python3 clauset_refit.py DATEI1 LABEL1 [DATEI2 LABEL2 ...] \
        [--outdir results/statistik/fit_gepoolt] [--bootstrap 500] [--gof 500] \
        [--fixed-xmin 4] [--jobs 8] [--seed 12345]
"""
import argparse, json, os, sys
import numpy as np

# powerlaw ist gespraechig; Warnungen der Fit-Routine unterdruecken
import warnings
warnings.filterwarnings("ignore")
np.seterr(all="ignore")
import powerlaw

# (Schluessel, Anzeigename, genestet?)
# "Potenzgesetz mit Cutoff" enthaelt das reine Potenzgesetz als Grenzfall
# (lambda -> 0) -> genesteter Test: p aus chi^2, R als rohe Log-Likelihood-
# Differenz. Die uebrigen Alternativen sind nicht genestet -> Vuong-Test mit
# normiertem R (Clauset Gl. C.6).
ALTERNATIVES = [
    ("lognormal",           "Lognormal",               False),
    ("exponential",         "Exponential",             False),
    ("stretched_exponential", "Stretched Exponential", False),
    ("truncated_power_law", "Potenzgesetz mit Cutoff", True),
]


def load(fn):
    """Liest eine Clustergroessen-CSV.

    Akzeptiert eine Spalte (`size`) ebenso wie zwei Spalten (`task,size`) --
    im zweiten Fall wird die letzte Spalte genommen.
    """
    d = np.atleast_1d(np.genfromtxt(fn, delimiter=",", skip_header=1))
    if d.ndim > 1:
        d = d[:, -1]
    d = d[np.isfinite(d)]
    d = d[d > 0]
    return np.sort(d.astype(int))


def fit_of(data, xmin=None):
    return powerlaw.Fit(data, discrete=True, xmin=xmin, verbose=False,
                        estimate_discrete=False)


def pl_stats(fit):
    """(S, x_min, sigma, KS-D, n_tail) -- robust auch bei fest vorgegebenem x_min
    (dann setzt powerlaw die Kurzform-Attribute des Fit-Objekts nicht)."""
    pl = fit.power_law
    return (float(pl.alpha), float(fit.xmin), float(pl.sigma),
            float(pl.D), int(pl.n))


# --------------------------------------------------------------------------
# Worker fuer die beiden Bootstraps (Modulebene -> picklebar fuer Multiprocessing)
# --------------------------------------------------------------------------
_BOOT = {}


def _boot_init(data, alpha, xmin, ntail, seed0):
    _BOOT["data"] = data
    _BOOT["alpha"] = alpha
    _BOOT["xmin"] = xmin
    _BOOT["ntail"] = ntail
    _BOOT["seed0"] = seed0


def _boot_param(i):
    """Nichtparametrischer Bootstrap: Resample der Daten, x_min neu schaetzen."""
    rng = np.random.default_rng(_BOOT["seed0"] + i)
    d = _BOOT["data"]
    s = rng.choice(d, size=len(d), replace=True)
    try:
        a_, x_, sg_, D_, nt_ = pl_stats(fit_of(s))
        return a_, x_, nt_, D_
    except Exception:
        return None


def _discrete_pl_sample(rng, n, alpha, xmin):
    """Diskrete Potenzgesetz-Zufallszahlen (Clauset Anhang D, exakte Methode)."""
    if n <= 0:
        return np.empty(0, dtype=int)
    # CDF-Tabelle bis dorthin, wo die Restmasse vernachlaessigbar ist
    hi = int(xmin)
    tail = 1.0
    xs = [int(xmin)]
    from scipy.special import zeta
    norm = zeta(alpha, xmin)
    while tail > 1e-9 and hi < 10**7:
        hi *= 2
        xs = np.arange(int(xmin), hi + 1)
        pmf = xs.astype(float) ** (-alpha) / norm
        tail = 1.0 - pmf.sum()
    cdf = np.cumsum(xs.astype(float) ** (-alpha) / norm)
    cdf /= cdf[-1]
    u = rng.random(n)
    return xs[np.searchsorted(cdf, u)]


def _boot_gof(i):
    """Semi-parametrischer Bootstrap fuer den p-Wert (Clauset Abschn. 4.1)."""
    rng = np.random.default_rng(_BOOT["seed0"] + 10_000_000 + i)
    d = _BOOT["data"]
    n = len(d)
    ptail = _BOOT["ntail"] / n
    below = d[d < _BOOT["xmin"]]
    k = int(rng.binomial(n, ptail))
    tail = _discrete_pl_sample(rng, k, _BOOT["alpha"], _BOOT["xmin"])
    if len(below) and n - k > 0:
        rest = rng.choice(below, size=n - k, replace=True)
    else:
        rest = np.empty(0, dtype=int)
    syn = np.concatenate([tail, rest]).astype(int)
    try:
        return pl_stats(fit_of(syn))[3]
    except Exception:
        return None


def mle_unbounded(data, xmin):
    """Diskrete MLE fuer alpha ohne Parametergrenze.

    powerlaw sucht alpha nur in [0, 3] und meldet bei steileren Verteilungen
    stumpf 3.000. Fuer den x_min-Scan ist das fatal, weil dort systematisch
    grosse Exponenten auftreten -- der Scan wuerde ein Plateau bei 3 zeigen,
    das es nicht gibt.
    """
    from scipy.optimize import minimize_scalar
    from scipy.special import zeta
    d = data[data >= xmin]
    n = len(d)
    ls = np.log(d).sum()
    r = minimize_scalar(lambda a: n * np.log(zeta(a, xmin)) + a * ls,
                        bounds=(1.001, 30.0), method="bounded",
                        options=dict(xatol=1e-6))
    return float(r.x)


def xmin_scan(data, n_points=40):
    """KS-Abstand D und Exponent S als Funktion von x_min.
    Zeigt, ob es ueberhaupt einen Bereich gibt, in dem das Potenzgesetz passt."""
    cands = np.unique(data)
    cands = cands[(data[:, None] >= cands).sum(axis=0) >= 30]
    if len(cands) > n_points:
        # logarithmisch abtasten: der Bereich kleiner x_min ist der interessante
        idx = np.unique(np.round(np.logspace(0, np.log10(len(cands)), n_points)
                                 ).astype(int) - 1)
        idx = idx[(idx >= 0) & (idx < len(cands))]
        cands = cands[idx]
    rows = []
    for xm in cands:
        try:
            a_, x_, sg_, D_, nt_ = pl_stats(fit_of(data, xmin=float(xm)))
            if abs(a_ - 3.0) < 1e-3:          # Paketgrenze erreicht -> frei nachrechnen
                a_ = mle_unbounded(data, float(xm))
            rows.append((float(xm), a_, D_, nt_))
        except Exception:
            continue
    return rows


def analyse(fn, label, args):
    data = load(fn)
    n = len(data)
    print(f"\n{'='*72}\n{label}  ({fn})\n{'='*72}")
    print(f"n = {n} Cluster, min = {data.min()}, max = {data.max()}, "
          f"Median = {np.median(data):.0f}")

    fit = fit_of(data)
    alpha, xmin, sigma, D, ntail = pl_stats(fit)
    print(f"\n-- Clauset-Fit (x_min via KS-Minimierung, MLE diskret) --")
    print(f"   S (alpha) = {alpha:.4f}   sigma_MLE = {sigma:.4f}")
    print(f"   x_min     = {xmin:.0f}    n_tail = {ntail} ({100*ntail/n:.1f} % der Daten)")
    print(f"   KS-Abstand D = {D:.4f}")

    # Naeherungsformel (Berichts-Gleichung, xmin - 1/2) als Kreuzcheck
    fit_approx = powerlaw.Fit(data, discrete=True, xmin=xmin, verbose=False,
                              estimate_discrete=True)
    print(f"   [Kreuzcheck Naeherung 1+n/sum ln(x/(x_min-1/2)): "
          f"S = {float(fit_approx.power_law.alpha):.4f}]")

    res = dict(label=label, file=fn, n=n, x_max=int(data.max()),
               alpha=alpha, sigma_mle=sigma, xmin=xmin, n_tail=ntail, D=D,
               alpha_approx=float(fit_approx.power_law.alpha))

    # ---- Bootstrap-CI ----------------------------------------------------
    from concurrent.futures import ProcessPoolExecutor
    if args.bootstrap > 0:
        print(f"\n-- Bootstrap ({args.bootstrap} Resamples, x_min jedes Mal neu) --")
        with ProcessPoolExecutor(max_workers=args.jobs, initializer=_boot_init,
                                 initargs=(data, alpha, xmin, ntail, args.seed)) as ex:
            out = [r for r in ex.map(_boot_param, range(args.bootstrap), chunksize=4)
                   if r is not None]
        ba = np.array([o[0] for o in out])
        bx = np.array([o[1] for o in out])
        lo, hi = np.percentile(ba, [2.5, 97.5])
        print(f"   S      = {alpha:.3f}  [{lo:.3f}, {hi:.3f}]  (95 %, Perzentil)")
        print(f"   SD(S)  = {ba.std(ddof=1):.4f}   (vgl. sigma_MLE = {sigma:.4f})")
        print(f"   x_min  Median = {np.median(bx):.0f}, "
              f"95 %-Bereich [{np.percentile(bx,2.5):.0f}, {np.percentile(bx,97.5):.0f}]")
        res.update(boot_n=len(ba), alpha_ci=[float(lo), float(hi)],
                   alpha_sd_boot=float(ba.std(ddof=1)),
                   xmin_median_boot=float(np.median(bx)),
                   xmin_ci=[float(np.percentile(bx, 2.5)), float(np.percentile(bx, 97.5))],
                   boot_alphas=ba.tolist())

    # ---- Goodness-of-Fit -------------------------------------------------
    if args.gof > 0:
        print(f"\n-- Goodness-of-Fit ({args.gof} synthetische Datensaetze) --")
        with ProcessPoolExecutor(max_workers=args.jobs, initializer=_boot_init,
                                 initargs=(data, alpha, xmin, ntail, args.seed)) as ex:
            Ds = [r for r in ex.map(_boot_gof, range(args.gof), chunksize=4)
                  if r is not None]
        Ds = np.array(Ds)
        p_gof = float((Ds >= D).mean())
        print(f"   D_beobachtet = {D:.4f},  Median D_synthetisch = {np.median(Ds):.4f}")
        print(f"   p = {p_gof:.3f}  -> {'Potenzgesetz plausibel' if p_gof > 0.1 else 'Potenzgesetz verworfen'}"
              f"  (Clauset-Kriterium: p > 0.1)")
        res.update(gof_n=len(Ds), p_gof=p_gof, D_synth_median=float(np.median(Ds)))

    # ---- Likelihood-Ratio-Tests -----------------------------------------
    print(f"\n-- Likelihood-Ratio-Tests (alle oberhalb x_min = {xmin:.0f}) --")
    print(f"   {'Alternative':<26} {'R':>10} {'p':>8}   Interpretation")
    comps = {}
    for key, name, nested in ALTERNATIVES:
        try:
            if nested:
                R, p = fit.distribution_compare("power_law", key, nested=True)
            else:
                R, p = fit.distribution_compare("power_law", key, nested=False,
                                                normalized_ratio=True)
        except Exception as e:
            print(f"   {name:<26} -- Fit fehlgeschlagen ({e})")
            continue
        if p > 0.1:
            verdict = "nicht unterscheidbar"
        elif R > 0:
            verdict = "Potenzgesetz bevorzugt"
        else:
            verdict = "Alternative bevorzugt"
        print(f"   {name:<26} {R:>10.3f} {p:>8.4f}   {verdict}"
              f"{'  [genestet, chi^2]' if nested else ''}")
        comps[key] = dict(name=name, R=float(R), p=float(p), verdict=verdict,
                          nested=bool(nested))
    res["compare"] = comps

    # Parameter der Alternativen (fuer die Diskussion)
    try:
        res["lognormal_mu"] = float(fit.lognormal.mu)
        res["lognormal_sigma"] = float(fit.lognormal.sigma)
        res["exponential_lambda"] = float(fit.exponential.Lambda)
        res["tpl_alpha"] = float(fit.truncated_power_law.alpha)
        res["tpl_lambda"] = float(fit.truncated_power_law.Lambda)
        lam = res["tpl_lambda"]
        cut = f"{1/lam:.0f}" if lam > 0 else "unendlich (lambda -> 0)"
        print(f"   [Cutoff-Fit: alpha = {res['tpl_alpha']:.3f}, "
              f"lambda = {lam:.3e} -> Cutoff bei n ~ {cut}]")
    except Exception:
        pass

    # ---- Kontrollfit bei festem x_min ------------------------------------
    if args.fixed_xmin:
        fx = float(args.fixed_xmin)
        fa, fxm, fsg, fD, fnt = pl_stats(fit_of(data, xmin=fx))
        print(f"\n-- Kontrolle: festes x_min = {fx:.0f} (Bericht) --")
        print(f"   S = {fa:.4f} +/- {fsg:.4f}   n_tail = {fnt}   D = {fD:.4f}")
        res["fixed"] = dict(xmin=fx, alpha=fa, sigma=fsg, n_tail=fnt, D=fD)
    res["xmin_scan"] = xmin_scan(data)
    res["_fit"] = fit
    res["_data"] = data
    return res


def make_plots(results, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ["#185FA5", "#D85A30", "#2E8B57"]

    # (1) CCDF mit Potenzgesetz + Alternativen
    fig, axes = plt.subplots(1, len(results), figsize=(6.2*len(results), 5.0),
                             squeeze=False)
    for k, r in enumerate(results):
        ax = axes[0][k]
        fit = r["_fit"]
        fit.plot_ccdf(ax=ax, color="0.35", marker=".", ms=4, linestyle="none",
                      label="Daten (CCDF)")
        fit.power_law.plot_ccdf(ax=ax, color=colors[0], lw=2.2,
                                label=f"Potenzgesetz  S={r['alpha']:.2f}")
        fit.lognormal.plot_ccdf(ax=ax, color=colors[1], lw=1.6, ls="--",
                                label="Lognormal")
        fit.truncated_power_law.plot_ccdf(ax=ax, color=colors[2], lw=1.6, ls="-.",
                                          label="Potenzgesetz + Cutoff")
        fit.exponential.plot_ccdf(ax=ax, color="#7B3FA0", lw=1.6, ls=":",
                                  label="Exponential")
        ax.axvline(r["xmin"], color="gray", ls=":", lw=1)
        ax.text(r["xmin"]*1.08, 1.5e-3, f"$x_{{min}}$={r['xmin']:.0f}",
                color="gray", fontsize=9)
        ci = r.get("alpha_ci")
        sub = (f"S = {r['alpha']:.2f} [{ci[0]:.2f}, {ci[1]:.2f}]" if ci
               else f"S = {r['alpha']:.2f}")
        pg = r.get("p_gof")
        if pg is not None:
            sub += f",  GoF p = {pg:.2f}"
        ax.set_title(f"{r['label']}\n{sub}", fontsize=11)
        ax.set_xlabel("Clustergroesse $n$")
        ax.set_ylabel(r"$P(N \geq n)$")
        ax.grid(alpha=.3, which="both")
        ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("Clauset-Refit: Potenzgesetz gegen Alternativen", fontsize=12)
    fig.tight_layout()
    p1 = os.path.join(outdir, "clauset_ccdf_alternatives.png")
    fig.savefig(p1, dpi=150); plt.close(fig)

    # (2) Bootstrap-Verteilung von S
    boots = [r for r in results if r.get("boot_alphas")]
    if boots:
        fig, ax = plt.subplots(figsize=(7.0, 4.6))
        for k, r in enumerate(boots):
            ba = np.array(r["boot_alphas"])
            ax.hist(ba, bins=40, alpha=.55, color=colors[k % len(colors)],
                    label=f"{r['label']}: {r['alpha']:.2f} "
                          f"[{r['alpha_ci'][0]:.2f}, {r['alpha_ci'][1]:.2f}]")
            ax.axvline(r["alpha"], color=colors[k % len(colors)], lw=2)
        ax.axvline(1.63, color="k", ls="--", lw=1.5,
                   label="Sand et al. 2013 (3D): S = 1.63")
        ax.set_xlabel("Potenzgesetz-Exponent $S$")
        ax.set_ylabel("Haeufigkeit (Bootstrap)")
        ax.set_title("Bootstrap-Verteilung des Exponenten (x_min je Resample neu)")
        ax.legend(fontsize=9); ax.grid(alpha=.3)
        fig.tight_layout()
        p2 = os.path.join(outdir, "clauset_bootstrap_alpha.png")
        fig.savefig(p2, dpi=150); plt.close(fig)
    else:
        p2 = None

    # (3) x_min-Diagnose: D(x_min) und S(x_min)
    scans = [r for r in results if r.get("xmin_scan")]
    if scans:
        fig, (axD, axS) = plt.subplots(1, 2, figsize=(12.0, 4.6))
        for k, r in enumerate(scans):
            xs = np.array([q[0] for q in r["xmin_scan"]])
            al = np.array([q[1] for q in r["xmin_scan"]])
            Ds = np.array([q[2] for q in r["xmin_scan"]])
            c = colors[k % len(colors)]
            axD.semilogx(xs, Ds, "o-", ms=4, color=c, label=r["label"])
            axS.semilogx(xs, al, "o-", ms=4, color=c, label=r["label"])
            axD.plot([r["xmin"]], [r["D"]], "*", ms=16, color=c)
        axD.set_xlabel("$x_{min}$"); axD.set_ylabel("KS-Abstand $D$")
        axD.set_title("KS-Abstand steigt mit $x_{min}$\n-> kein Bereich, in dem das Potenzgesetz besser passt")
        axS.axhline(1.63, color="k", ls="--", lw=1.2, label="Sand et al. (3D): 1.63")
        axS.set_xlabel("$x_{min}$"); axS.set_ylabel("Exponent $S$")
        axS.set_title("Exponent haengt stark von $x_{min}$ ab\n(Stern = KS-Optimum)")
        for ax in (axD, axS):
            ax.grid(alpha=.3, which="both"); ax.legend(fontsize=9)
        fig.suptitle("$x_{min}$-Diagnose (jeder Punkt = eigener MLE-Fit)", fontsize=12)
        fig.tight_layout()
        p3 = os.path.join(outdir, "clauset_xmin_scan.png")
        fig.savefig(p3, dpi=150); plt.close(fig)
    else:
        p3 = None
    return p1, p2, p3


def write_report(results, outdir, args):
    md = os.path.join(outdir, "clauset_refit.md")
    with open(md, "w") as f:
        f.write("# Clauset-Refit der Kaskaden-Clustergroessen\n\n")
        f.write(f"Paket: `powerlaw` {powerlaw.__version__} (Alstott et al. 2014), "
                f"Methodik: Clauset, Shalizi & Newman (2009).\n")
        f.write(f"Bootstrap: {args.bootstrap} Resamples, "
                f"Goodness-of-Fit: {args.gof} synthetische Datensaetze, Seed {args.seed}.\n\n")

        f.write("## Hauptergebnis\n\n")
        f.write("| Datensatz | n | x_min | n_tail | S | sigma_MLE | 95 %-CI (Bootstrap) | KS D | GoF p |\n")
        f.write("|---|---:|---:|---:|---:|---:|---|---:|---:|\n")
        for r in results:
            ci = r.get("alpha_ci")
            cis = f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "--"
            pg = f"{r['p_gof']:.3f}" if r.get("p_gof") is not None else "--"
            f.write(f"| {r['label']} | {r['n']} | {r['xmin']:.0f} | {r['n_tail']} | "
                    f"{r['alpha']:.3f} | {r['sigma_mle']:.3f} | {cis} | {r['D']:.4f} | {pg} |\n")

        f.write("\n## Likelihood-Ratio-Tests (Potenzgesetz vs. Alternative)\n\n")
        f.write("R > 0 = Potenzgesetz bevorzugt, R < 0 = Alternative bevorzugt; "
                "p > 0.1 = Unterschied statistisch nicht signifikant.\n"
                "Nicht genestete Alternativen: Vuong-Test mit normiertem R. "
                "Genestet (Cutoff): R = rohe Log-Likelihood-Differenz, p aus chi^2.\n\n")
        f.write("| Datensatz | Alternative | Test | R | p | Interpretation |\n")
        f.write("|---|---|---|---:|---:|---|\n")
        for r in results:
            for key, c in r.get("compare", {}).items():
                f.write(f"| {r['label']} | {c['name']} | "
                        f"{'genestet' if c.get('nested') else 'Vuong'} | {c['R']:.3f} | "
                        f"{c['p']:.4f} | {c['verdict']} |\n")

        if any("fixed" in r for r in results):
            f.write(f"\n## Kontrollfit bei festem x_min = {args.fixed_xmin} "
                    f"(Wert aus dem Projektbericht)\n\n")
            f.write("| Datensatz | S | sigma | n_tail | KS D |\n|---|---:|---:|---:|---:|\n")
            for r in results:
                fx = r.get("fixed")
                if fx:
                    f.write(f"| {r['label']} | {fx['alpha']:.3f} | {fx['sigma']:.3f} | "
                            f"{fx['n_tail']} | {fx['D']:.4f} |\n")

        f.write("\n## x_min-Diagnose (KS-Abstand als Funktion von x_min)\n\n")
        for r in results:
            sc = r.get("xmin_scan")
            if not sc:
                continue
            f.write(f"\n**{r['label']}**\n\n| x_min | n_tail | S | KS D |\n|---:|---:|---:|---:|\n")
            for xm, al, D_, nt in sc[:14]:
                f.write(f"| {xm:.0f} | {nt} | {al:.3f} | {D_:.4f} |\n")

        f.write("\n## Foliensatz-Formulierung\n\n")
        for r in results:
            ci = r.get("alpha_ci")
            line = f"**{r['label']}:** S = {r['alpha']:.2f} +/- {r['sigma_mle']:.2f}"
            if ci:
                line += f" (95 %-Bootstrap-CI [{ci[0]:.2f}, {ci[1]:.2f}])"
            line += f", x_min = {r['xmin']:.0f} per KS-Minimierung"
            if r.get("p_gof") is not None:
                line += f", Goodness-of-Fit p = {r['p_gof']:.2f}"
            f.write(f"- {line}\n")

    # JSON ohne die nicht serialisierbaren Objekte
    js = []
    for r in results:
        d = {k: v for k, v in r.items() if not k.startswith("_")}
        js.append(d)
    with open(os.path.join(outdir, "clauset_refit.json"), "w") as f:
        json.dump(js, f, indent=2)
    return md


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs", nargs="+", help="DATEI LABEL [DATEI LABEL ...]")
    ap.add_argument("--outdir", default="results/statistik/fit_gepoolt")
    ap.add_argument("--bootstrap", type=int, default=500)
    ap.add_argument("--gof", type=int, default=500)
    ap.add_argument("--fixed-xmin", type=float, default=4.0)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--seed", type=int, default=12345)
    a = ap.parse_args()
    if len(a.pairs) % 2:
        sys.exit("Bitte Paare aus DATEI und LABEL angeben.")
    os.makedirs(a.outdir, exist_ok=True)

    results = []
    for i in range(0, len(a.pairs), 2):
        results.append(analyse(a.pairs[i], a.pairs[i+1], a))

    p1, p2, p3 = make_plots(results, a.outdir)
    md = write_report(results, a.outdir, a)
    print(f"\ngeschrieben:\n  {md}\n  {os.path.join(a.outdir,'clauset_refit.json')}\n  {p1}")
    for extra in (p2, p3):
        if extra:
            print(f"  {extra}")


if __name__ == "__main__":
    main()
