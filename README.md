# Was it the Sun?

A reproducible rebuild of a climate figure circulating on X which argues that the
warming since the Little Ice Age is solar in origin.

The rebuild keeps the original's visual architecture — one wide panel, 1600–2025,
solar-minimum shading, annotated turning points — so the two can be read side by
side. What changes is the parts that were wrong.

![The corrected figure](figure_climate_1600_2025.png)

## What the original gets wrong

**1. Raw TSI is plotted as if it were a radiative forcing.** Total solar irradiance
in W/m² at the top of the atmosphere is not a forcing. The Sun illuminates a
cross-section πr² while the Earth radiates from 4πr², and about 30 % is reflected,
so the conversion is `ΔTSI/4 × (1−α)` = **0.175**, and skipping it overstates the
solar term by its reciprocal, **5.7×**, on geometry alone. Both numbers follow from
`ALBEDO` in `climate_figure.py` and are computed there rather than typed into the
figure; the two quoted here are checked against that constant by `test_pipeline.py`.
Scaling the solar axis to the 11-year cycle peaks compounds the error.

**2. HadCET is presented as a global record.** The sharp warming around 1730 is a
Central England signal. Measured identically — mean of 1680–1700 against mean of
1725–1735 — England warms **+1.21 K** while the Northern Hemisphere warms
**+0.049 K**. The English excursion is the **98th percentile** of all 311 55-year
windows in the CET record; the hemispheric one barely moves.

A "55-year window" starting at year *y0* means the mean over *y0…y0+20* subtracted
from the mean over *y0+45…y0+55*. The published comparison is the *y0* = 1680
member of that family. The definition lives in one place in `climate_figure.py`, and
`audit.py` calls the same function, so the figure and its audit count the same 311
windows (start years 1659–1969) rather than two slightly different populations.

The *ratio* of those two numbers is 24.5, but this repository does not lead with
it. Its denominator is a hemispheric mean close to zero, which makes the quotient
unstable — across 54 sliding windows the same ratio runs from 0.1 to 43.2, and 6 %
of them exceed our value without anything remarkable happening. Two absolute
numbers and a percentile are the defensible statement; a large quotient is not.
(`audit.py` block 8 computes these.)

**3. The Spörer Minimum (~1460–1550) is labelled on an axis that starts in 1600,**
roughly 150 years from the period it names.

## What this repository finds

Every value annotated in the figure is computed by the pipeline from the validated
files in the same run. The rounded values below are transcribed into this README by
hand and reproduce the current pinned run; `test_pipeline.py` checks the headline
figures against `derive()` so they cannot drift out of step unnoticed.

Both forcings are read at **2018**, the last year both are defined: the CO₂ record
runs to 2025, but a *centred* 11-year solar mean necessarily stops five years short
of the end of the TSI record. Comparing CO₂ at 2025 with a solar value at 2018 and
calling both "today" would inflate every ratio below.

| | |
|---|---|
| CO₂ forcing, 1750→2018 | **+2.085 W/m²** |
| Solar forcing, 1750→2018 | **+0.108 W/m²** |
| Ratio | **19 : 1** |
| Ratio, measured from the Maunder trough (1692→2018) | **12 : 1** |
| Solar forcing trend since 1980 | **−0.009 W/m²/decade** (to 2018) |
| HadCRUT5 trend since 1980 | **+0.211 K/decade** |
| UAH lower-troposphere trend since 1980 | **+0.157 K/decade** |
| Warming 1680–1700 → 1725–1735, HadCET | **+1.206 K** (single series) |
| …same window, N. Hemisphere reconstruction | **+0.049 K**, 5–95 % **+0.015 … +0.093**, 99 % of members positive |
| …same window, S. Hemisphere reconstruction | **+0.035 K**, 5–95 % **−0.103 … +0.209**, only 69 % positive |
| HadCET excursion vs all 311 CET 55-year windows (1659–1969) | **98th percentile** |

The hemispheric figures are medians over 100 ensemble members, with the 5–95 %
spread of the same members. The Southern Hemisphere interval spans zero, so **the
sign of the SH excursion is not established** and is not claimed here. The Northern
Hemisphere result is what carries the regional-versus-global argument.

Two results carry most of the weight:

**The Sun and the Earth have moved in opposite directions since 1980** — within
NRLTSI2, the only TSI dataset loaded here. No choice of *axis* can repair a
disagreement in sign, and this result is far less sensitive to the reconstruction
than the amplitude comparison is; it is not, however, independent of it. See
[the objection this repository cannot settle](#the-objection-this-repository-cannot-settle).
The strongest 11-year mean TSI in the record is **1955**, not today; 2018 ranks 35th
of 404 years.

**The source figure's own reference (3) does not support it.** Thompson et al.
(1986), the Quelccaya ice cap in Peru, is cited there for "seen in many locations
globally". Over the same 1680–1700 → 1725–1735 window its δ¹⁸O changes by −0.80 ‰,
and Central England correlates with it at r = −0.018 (annual, n = 142) over
1659–1800 — `audit.py` block 8 computes this alongside the smoothed value. The
paper's own abstract highlights 1590–1630, 1800–1840 and 1915–1940 — not 1690–1730.

Note what this does *not* say. Tropical δ¹⁸O is precipitation-influenced, not a
thermometer — the authors' companion title is "A 1500-Year Record of Tropical
Precipitation" — so the sign of that change is **not** read here as cooling, and no
proxy-to-temperature model is applied. The defensible statement is the narrower one:
the record provides no clear evidence of a parallel global warming, and therefore
cannot carry the claim made of it.

## Reproducing it

```bash
pip install -r requirements.txt
python fetch_data.py
python climate_figure.py
```

Exit code 0 with `figure_climate_1600_2025.png` written means every series passed
validation. Exit code 1 means one failed — read stderr and fix the data or the
config; do not weaken the check.

```bash
python robustness.py
python test_pipeline.py
```

`robustness.py` prints the sensitivity of the CO₂:solar ratio to albedo, to the
reference year, and to skipping the geometric conversion altogether, plus a
cross-check of the CO₂ forcing against IPCC AR6. `test_pipeline.py` asserts the
guarantees the figure depends on: that both forcings are read at the same year, that
the ensemble intervals are recomputed from the ensemble files, that the plotted bar
and its whiskers are the same statistic, that every MANIFEST field is present and
non-empty, and that no series is loaded without being used.

## What this repository claims, and what it does not

**It claims** that one specific figure does not support the conclusion printed on
it: the solar term it draws is not a radiative forcing, the 1730 excursion it
presents as global is a Central England signal, and the reference it cites for
global synchrony does not show what it is cited for.

Each series used here is loaded from its explicitly named published archive. Where
the source figure cannot be reproduced directly — most notably GloSAT, which needs a
free CEDA account — the replacement or additional comparison series is identified
under its real name and is never presented as the original dataset. HadCRUT5 is
labelled HadCRUT5, not GloSAT; the Neukom 2018 hemispheric reconstructions are an
additional comparison the source figure does not use at all.

**It does not claim** to settle attribution of modern warming. That question is
answered by a literature this repository does not attempt to summarise — fingerprint
studies, energy-balance constraints, ocean heat content, model–observation
comparison. Nothing here depends on that literature, and nothing here replaces it.

**It does not claim** its own numbers are the last word. The window definition, the
solar reconstruction and the pre-industrial CO₂ reference are all contestable
choices; each is quantified below rather than defended, and `audit.py` exists to
attack them. Where a choice turned out to be doing the work — the excursion ratio —
the figure was changed rather than the wording.

## What was checked against itself

`audit.py` attacks the conclusions above with the loaded data. The checks it runs,
and what each one did to the figure:

- **CO₂ splice** — the step across the 1959 ice-core → Mauna Loa join is +0.78 ppm,
  smaller than the +0.88 ppm/yr rise of the 1960s. The splice does not manufacture a
  jump. Mauna Loa runs ~2.3 ppm above the Law Dome spline over their overlap, as
  expected for a Hawaiian station against a smoothed Antarctic core.
- **TSI time axis** — the "days since 1610-01-01" conversion puts cycle maxima at
  1957, 1981, 2000 and 2015 against expected ~1958, ~1979, ~2000–2002 and ~2014.
  Scatter of ±2 years in both directions, no systematic offset.
- **Reconstruction damping** — over the 1850–1999 overlap the NH reconstruction
  carries 93 % of HadCRUT5's standard deviation, so there is no sign of strong
  general variance compression. This check cannot rule out period-specific or
  frequency-dependent damping in the early 18th century, where no instrumental
  record exists to compare against. The ensemble spread, not this ratio, remains the
  operative statement of historical uncertainty.
- **CO₂ reference value** — the Law Dome spline gives 276.80 ppm at 1750 against the
  IPCC AR6 reference of 278.3 ppm. Using the published value instead would *lower*
  our CO₂ forcing by 0.029 W/m² and the ratio from 19.2 to about 18.9. The loaded
  value is kept because typing in a published constant is exactly what this pipeline
  forbids, and the difference is disclosed rather than pocketed.
- **Window choice** — across 45 alternative definitions (start year 1670–1690,
  baseline 15/20/25 years, peak offset +40–50/+45–55/+50–60) the hemispheric
  excursion never exceeds +0.12 K, median +0.007 K. The English one runs +0.18 to
  +1.26 K, median +0.84 K. **The published window is at the favourable end for
  England** — that is stated on the figure itself rather than left to be found. The
  hemispheric result, which is what the argument rests on, is insensitive to the
  choice.
- **Numbers quoted in this README** — block 8 recomputes the two figures that had
  no live source (the Quelccaya–HadCET correlation and the sliding-window ratio
  distribution), so nothing quoted here rests on a script that no longer exists.
- **The excursion ratio did not survive** and the figure was changed. See above.

### The objection this repository cannot settle

The solar reconstruction used here, NRLTSI2, is a low-variability one, and that
choice is contestable — higher-amplitude reconstructions of the Maunder-to-present
irradiance change exist. This repository does not adjudicate that debate. It
computes what the objection would have to be worth: for the solar forcing to *equal*
the CO₂ forcing, ΔTSI from 1750 to 2018 would have to be about **19× larger than
NRLTSI2 gives — near 11.9 W/m², roughly 1 % of the solar constant**. A critic who
rejects NRLTSI2 needs a reconstruction that reaches that, not merely a higher one.

The post-1980 sign result is **less** dependent on that choice, but it is not
independent of it. NRLTSI2 is calibrated against direct space-based observations
over the satellite era and reproduces their variability, so its modern trend is far
better constrained than its Maunder-era amplitude. It remains a proxy model built on
sunspot and facular indices regressed onto those observations, not a satellite
composite. This pipeline loads NRLTSI2 and nothing else, so what it can honestly say
is narrower than "reconstruction-independent": **within the NRLTSI2 dataset the trend
since 1980 is negative**, and no alternative TSI composite has been loaded here to
test that against.

## How it is kept honest

`data/` is not committed; `fetch_data.py` downloads the exact files in the MANIFEST.
Every series carries a MANIFEST entry with citation, URL, access date and the
**SHA-256 of the file actually used** — that hash has one home, and `EXPECT` does
not keep a second copy of it, because two hand-maintained copies drift and a
drifted hash is worse than none: it still looks checked. Validation and the
downloader both read `MANIFEST[key]["sha256"]`.

Checked before anything is plotted:

- **sha256** from the MANIFEST — a wrong file spanning the right years passes a year
  check in silence
- **coverage** — first and last year (`EXPECT`)
- **value_range** — catches unit errors and column mix-ups (`EXPECT`)
- **n_min** — catches a file that parsed into mostly NaN (`EXPECT`)
- **column names** — asserted against the real file header, so `CFG` cannot drift
  into documentation

MANIFEST completeness is enforced in code: an undocumented series or an orphan entry
raises `DataError`. A QC table prints for every series before any figure is written.

The rules the pipeline is held to are in [AGENTS.md](AGENTS.md), and they apply to
this figure too. Notably: if this figure drew the solar term smaller than its
computed forcing warrants, that would be the same error as the original with the
sign reversed. The solar series therefore also appears at its own scale, with the
magnification stated, so its real structure — the Maunder dip, Dalton, the 1955 peak
— stays visible.

## Sources

HadCET (Manley 1974; Parker et al. 1992) · HadCRUT5 (Morice et al. 2021) ·
UAH v6.1 (Spencer, Christy & Braswell 2017) · Law Dome CO₂ (MacFarling Meure et al.
2006) · NOAA GML Mauna Loa · NRLTSI2 (Coddington et al. 2016) ·
Quelccaya (Thompson et al. 1986) · Neukom et al. 2018 hemispheric reconstructions
from PAGES 2k v2.0.0 proxies.

Full citations, URLs, access dates and file hashes: the MANIFEST block at the top of
[climate_figure.py](climate_figure.py).

The SILSO sunspot number was loaded in an earlier revision but never plotted or used in any computation. It has been removed rather than left half-registered: it is a proxy for the same activity NRLTSI2 already reconstructs, and it cannot be converted to a radiative forcing. The reasoning is recorded in the MANIFEST block.

GloSAT is deliberately not used — it requires a free CEDA account. The global series
here is HadCRUT5 and is labelled HadCRUT5.

## Licence

Code, figure and documentation in this repository: MIT, see [LICENSE](LICENSE).

The **data is not covered by that licence and is not redistributed here**. Each
source keeps its own terms and its own citation requirements; `fetch_data.py`
downloads the files directly from the providers listed in the MANIFEST. If you use
the figure or reproduce the analysis, cite the underlying datasets, not this
repository — the citations are in the MANIFEST block of `climate_figure.py`.

The figure that prompted this rebuild is not reproduced here. It is identified in
the source header of `climate_figure.py` and its claims are quoted as text so they
can be checked against the data.
