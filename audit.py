"""Adversarial checks against this pipeline's own conclusions.

Each block below is an attempt to break a claim the figure makes, using only the
loaded data. Anything that survives is worth keeping; anything that does not has
to change the figure, not the wording.
"""

import numpy as np
import pandas as pd

import climate_figure as cf

series, _, _ = cf.load_all()
d = cf.derive(series)
W = "=" * 78


def head(n, s):
    print(f"\n{W}\n{n}. {s}\n{W}")


# --------------------------------------------------------------------------
head(1, "CO2 SPLICE: is the ice core / Mauna Loa join a real discontinuity?")
ice, mlo = series["co2_ice"], series["co2_mlo"]
ov = ice.index.intersection(mlo.index)
print(f"   overlap {int(ov.min())}-{int(ov.max())}, n={len(ov)}")
diff = (mlo.loc[ov] - ice.loc[ov])
print(f"   Mauna Loa minus Law Dome spline: mean {diff.mean():+.3f} ppm, "
      f"min {diff.min():+.3f}, max {diff.max():+.3f}")
sy = d["splice_year"]
print(f"   step across the splice at {sy}: "
      f"{float(mlo.loc[sy]) - float(ice.loc[sy - 1]):+.3f} ppm "
      f"(vs typical annual rise "
      f"{float(mlo.diff().loc[1960:1970].mean()):+.3f} ppm/yr in the 1960s)")
print("   VERDICT: a step much larger than the annual rise would mean the splice "
      "is\n   manufacturing a jump. Read the numbers above before trusting the "
      "CO2 curve.")


# --------------------------------------------------------------------------
head(2, "TSI TIME AXIS: does the converted year land on known solar maxima?")
raw = series["tsi"]
# Known sunspot-cycle maxima, used only as a ruler for the time axis.
for lo, hi, name in [(1955, 1962, "cycle 19 (max ~1958)"),
                     (1976, 1983, "cycle 21 (max ~1979)"),
                     (1998, 2004, "cycle 23 (max ~2000-2002)"),
                     (2011, 2016, "cycle 24 (max ~2014)")]:
    w = raw.loc[lo:hi]
    print(f"   {name:<28s} local TSI max at {int(w.idxmax())}")
print("   VERDICT: an off-by-one or a scaling error in 'days since 1610-01-01'\n"
      "   would show up as maxima landing consistently early or late.")


# --------------------------------------------------------------------------
head(3, "RECONSTRUCTION VARIANCE: is the hemispheric series damped?")
# Proxy reconstructions regress towards the mean; if the NH series has far less
# variance than an instrumental record over the same years, then a small NH
# excursion is partly an artefact of the method, not only of spatial averaging.
nh, had = d["nh"], d["hadcrut5"]
common = nh.index.intersection(had.index)
for lo, hi in [(1850, 1999), (1850, 1900)]:
    w = [x for x in common if lo <= x <= hi]
    a, b = nh.loc[w], had.loc[w]
    print(f"   {lo}-{hi}: sd(NH recon)={a.std():.4f} K   sd(HadCRUT5)={b.std():.4f} K"
          f"   ratio={a.std() / b.std():.2f}   r={np.corrcoef(a, b)[0, 1]:+.3f}")
print("   VERDICT: a ratio well below 1 would mean the reconstruction understates\n"
      "   real variability, making +0.05 K a lower bound. It does not, over this\n"
      "   overlap. But the overlap is 1850-1999 and the claim is about 1680-1735,\n"
      "   where no instrumental record exists to compare against: this test cannot\n"
      "   exclude period-specific or frequency-dependent damping there. The\n"
      "   ensemble spread remains the operative uncertainty statement.")


# --------------------------------------------------------------------------
head(4, "IS 24x REMARKABLE? Any single region beats a hemispheric mean.")
# Averaging kills variance. Part of the CET/NH ratio is arithmetic, not evidence.
cet, nhs = d["cet"], d["nh"]
common = cet.index.intersection(nhs.index)
w = [x for x in common if 1659 <= x <= 1900]
print(f"   1659-1900  sd(HadCET)={cet.loc[w].std():.4f} K  "
      f"sd(NH recon)={nhs.loc[w].std():.4f} K  ratio={cet.loc[w].std() / nhs.loc[w].std():.1f}")
# same statistic as the headline, but for every 55-year window, to see whether
# 1680-1735 is unusual or just typical of England vs a hemispheric mean
ratios = []
for y0 in range(1660, 1930, 5):
    a0, a1, b0, b1 = y0, y0 + 20, y0 + 45, y0 + 55
    if b1 > min(cet.index.max(), nhs.index.max()):
        break
    ce = cet.loc[b0:b1].mean() - cet.loc[a0:a1].mean()
    ne = nhs.loc[b0:b1].mean() - nhs.loc[a0:a1].mean()
    if abs(ne) > 1e-6:
        ratios.append(abs(ce / ne))
ratios = np.array(ratios)
print(f"   |CET excursion / NH excursion| over {len(ratios)} sliding windows: "
      f"median {np.median(ratios):.1f}, "
      f"5-95% [{np.percentile(ratios, 5):.1f}, {np.percentile(ratios, 95):.1f}]")
print(f"   our 1680-1735 window: {d['exc_ratio']:.1f}")
print("   VERDICT: if 24x sits near the middle of that distribution, the ratio is\n"
      "   ordinary regional-vs-hemispheric behaviour and must NOT be presented as\n"
      "   if it were specific to 1730.")


# --------------------------------------------------------------------------
head(5, "THE ABSOLUTE CLAIM: is the CET excursion itself unusual?")
# The defensible version of the argument is about the SIZE of the English event,
# not about a ratio that averaging guarantees.
# Uses the pipeline's own window function, so this audit cannot be measuring a
# different population from the one the published figure quotes.
ex = cf.excursion_windows(cet)
our = d["exc"]["cet"]
print(f"   CET {cf.EXC_LEN}-year excursions, {len(ex)} windows "
      f"(start years {int(ex.index.min())}-{int(ex.index.max())}): "
      f"median {ex.median():+.3f} K, 95th pct {ex.quantile(0.95):+.3f} K, "
      f"max {ex.max():+.3f} K at {int(ex.idxmax())}")
print(f"   our {cf.EXC_START}->{cf.EXC_START + cf.EXC_LEN} window: {our:+.3f} K  "
      f"= {(ex.to_numpy() < our).mean() * 100:.4f}th percentile")
print(f"   pipeline reports {d['cet_windows_n']} windows and "
      f"{d['cet_pct']:.4f}th percentile -- these must be identical")
assert len(ex) == d["cet_windows_n"], "audit and pipeline disagree on window count"
assert abs((ex.to_numpy() < our).mean() * 100 - d["cet_pct"]) < 1e-9


# --------------------------------------------------------------------------
head(6, "CO2 REFERENCE VALUE: does C0 match the published pre-industrial?")
c0 = float(series["co2_ice"].loc[1750])
print(f"   Law Dome spline at 1750: {c0:.2f} ppm")
print(f"   IPCC AR6 pre-industrial (1750) reference: 278.3 ppm")
print(f"   forcing difference from using 278.3 instead: "
      f"{cf.CO2_FORCING_COEFF * np.log(float(d['co2'].loc[d['compare_year']]) / 278.3) - d['f_co2_now']:+.4f} W/m2")
print("   VERDICT: a lower C0 inflates our CO2 forcing. If the gap matters at the\n"
      "   second decimal, say so rather than quietly keeping the friendlier value.")
print(f"\n{W}\n")


# --------------------------------------------------------------------------
head(7, "WINDOW CHOICE: does the result depend on 1680-1700 vs 1725-1735?")
# The published windows came from the brief, not from a source. Rather than
# defend the choice, vary it and see whether anything survives. If the numbers
# move a lot, the choice is doing the work and the claim has to be weakened.
cet_s, nh_s = d["cet"], d["nh"]


def exc_span(s, y0, base_len, peak_lo, peak_hi):
    a = s.loc[y0:y0 + base_len].mean()
    b = s.loc[y0 + peak_lo:y0 + peak_hi].mean()
    return float(b - a)


rows = []
for y0 in range(1670, 1691, 5):
    for base_len in (15, 20, 25):
        for peak_lo, peak_hi in ((40, 50), (45, 55), (50, 60)):
            rows.append((exc_span(cet_s, y0, base_len, peak_lo, peak_hi),
                         exc_span(nh_s, y0, base_len, peak_lo, peak_hi)))
ce = np.array([r[0] for r in rows])
ne = np.array([r[1] for r in rows])
print(f"   {len(rows)} window definitions (start 1670-1690, baseline 15/20/25 yr,")
print(f"   peak +40-50 / +45-55 / +50-60):")
print(f"     HadCET excursion: median {np.median(ce):+.3f} K, "
      f"range [{ce.min():+.3f}, {ce.max():+.3f}] K")
print(f"     NH   excursion:   median {np.median(ne):+.3f} K, "
      f"range [{ne.min():+.3f}, {ne.max():+.3f}] K")
print(f"     published choice: HadCET {d['exc']['cet']:+.3f} K, "
      f"NH {d['exc']['nh']:+.3f} K")
print(f"     HadCET stays above +0.5 K in {(ce > 0.5).mean() * 100:.0f}% of them; "
      f"NH stays below +0.2 K in {(ne < 0.2).mean() * 100:.0f}%.")
print("   VERDICT: if the English excursion collapses under other windows, the")
print("   published pair was doing the work and the claim must be weakened.")
print(f"\n{W}\n")


# --------------------------------------------------------------------------
head(8, "NUMBERS QUOTED IN THE README THAT NOTHING ELSE COMPUTES")
# Every figure in the README has to trace to a run. These two were transcribed
# from throwaway scripts and had no live source until now.
quel = series["quelccaya"]
both = pd.concat([quel.rename("q"), d["cet"].rename("cet")], axis=1).dropna()
for lo, hi in [(1659, 1800), (1659, 1984)]:
    w = both.loc[lo:hi]
    r_ann = float(np.corrcoef(w["q"], w["cet"])[0, 1])
    sm = w.rolling(11, center=True).mean().dropna()
    r_sm = float(np.corrcoef(sm["q"], sm["cet"])[0, 1])
    print(f"   Quelccaya d18O vs HadCET {lo}-{hi}: n={len(w)}  "
          f"r(annual)={r_ann:+.3f}  r(11-yr)={r_sm:+.3f}")
print("   (a correlation, not a temperature claim: see the MANIFEST note)")

# how ordinary is our excursion ratio among the sliding windows of block 4?
above = float((ratios > d["exc_ratio"]).mean() * 100)
print(f"\n   |CET/NH| ratio of the published window: {d['exc_ratio']:.1f}")
print(f"   share of the {len(ratios)} sliding windows that exceed it: {above:.0f}%")
print(f"   full range of those ratios: {ratios.min():.1f} to {ratios.max():.1f}")
print("   VERDICT: a value the ordinary distribution reaches this often is not")
print("   evidence on its own. This is why the ratio is not the headline.")
print(f"\n{W}\n")
