"""
Solar activity vs CO2 as drivers of observed warming, 1600-2025.

A corrected rebuild of a widely circulated figure (see Source.jpg) which argues
that the warming since the Little Ice Age is solar in origin. This script keeps
that figure's visual architecture -- one wide panel, 1600-2025, era shading,
annotated turning points -- so the two can be read side by side, and corrects
three defects in it:

  1. It plots HadCET (Central England) alongside global records on one
     temperature axis, so a regional excursion around 1730 reads as global.
  2. It labels the Spoerer Minimum (~1460-1550) at the left edge of an axis
     that starts in 1600, ~150 years from the period named.
  3. It plots raw TSI in W/m2 against temperature and CO2 on an axis scaled to
     the 11-year cycle peaks. Raw TSI is not a radiative forcing; the
     conversion is dTSI/4*(1-albedo) ~ 0.175. The effect is to overstate the
     solar term by roughly an order of magnitude.

Every value plotted here is loaded from a file in data/. Nothing is
reconstructed, and no series is smoothed except where a published smoothing is
used and labelled as such.

==============================================================================
MANIFEST
==============================================================================

hadcet
  Mean Central England Temperature, monthly and annual, 1659-present.
  Manley (1974) Q.J.R.Meteorol.Soc. 100:389-405; Parker et al. (1992)
  Int.J.Climatol. 12:317-342; Parker & Horton (2005) Int.J.Climatol. 25:1173.
  URL:      https://www.metoffice.gov.uk/hadobs/hadcet/
  File:     https://www.metoffice.gov.uk/hadobs/hadcet/data/meantemp_monthly_totals.txt
  Accessed: 2026-08-30
  Note:     Absolute degC. Regional (English Midlands), NOT a global record.
            Converted here to an anomaly on its own 1961-1990 mean.

hadcrut5
  HadCRUT5 analysis, global mean surface temperature anomaly, annual, 1850-.
  Morice et al. (2021) JGR Atmospheres 126:e2019JD032361.
  URL:      https://www.metoffice.gov.uk/hadobs/hadcrut5/
  File:     https://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.0.2.0/
            analysis/diagnostics/HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv
  Accessed: 2026-08-30
  Note:     Already on the 1961-1990 baseline. No re-referencing applied.
            The source figure labels its global series GloSAT; GloSAT requires
            a free CEDA account and is not used here. This series is labelled
            HadCRUT5 because that is what it is.

uah_lt
  UAH satellite lower-troposphere temperature anomaly v6.1, monthly, 1978-12-.
  Spencer, Christy & Braswell (2017) Asia-Pac.J.Atmos.Sci. 53:121-130.
  URL:      https://www.nsstc.uah.edu/climate/
  File:     https://www.nsstc.uah.edu/data/msu/v6.1/tlt/uahncdc_lt_6.1.txt
  Accessed: 2026-08-30
  Note:     Published on a 1991-2020 baseline. It CANNOT be re-referenced to
            1961-1990 because the record begins in 1978 -- there is no overlap
            with the baseline period. It is instead offset onto HadCRUT5 over
            their common years, and that offset is printed in the QC table.
            Lower troposphere is not a surface record; it is drawn dashed and
            labelled as a different quantity.

co2_ice
  Law Dome ice core CO2, published smoothing spline, 1-2004 AD.
  MacFarling Meure et al. (2006) Geophys.Res.Lett. 33:L14810;
  Etheridge et al. (1996) JGR 101:4115-4128.
  URL:      https://www.ncei.noaa.gov/access/paleo-search/study/9959
  File:     https://www.ncei.onaa.gov/pub/data/paleo/icecore/antarctica/law/law2006.txt
  Accessed: 2026-08-30
  Note:     The spline is the authors' published product (Enting 1987 method),
            not a smoothing applied here. Used only before the Mauna Loa record
            begins; the splice is drawn with both series visible in overlap.

co2_mlo
  NOAA GML Mauna Loa annual mean CO2, 1959-present.
  Keeling et al. (1976); NOAA GML, Lan, Tans & Thoning.
  URL:      https://gml.noaa.gov/ccgg/trends/
  File:     https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_mlo.txt
  Accessed: 2026-08-30

tsi
  NRLTSI2 total solar irradiance reconstruction, annual, 1610-present.
  Coddington et al. (2016) Bull.Amer.Meteor.Soc. 97:1265-1282;
  Lean (2000) Geophys.Res.Lett. 27:2425-2428.
  URL:      https://lasp.colorado.edu/lisird/data/nrl2_tsi_P1Y
  File:     https://lasp.colorado.edu/lisird/latis/dap/nrl2_tsi_P1Y.csv
  Accessed: 2026-08-30
  Note:     Time axis is days since 1610-01-01, converted to calendar year in
            the loader. Absolute scale ~1361 W/m2.

quelccaya
  Quelccaya Ice Cap summit core, annual d18O, 744-1984 AD. Peru, 13.9S, 5670 m.
  Thompson, Mosley-Thompson, Dansgaard & Grootes (1986) Science 234:361-364,
  "The Little Ice Age as Recorded in the Stratigraphy of the Tropical Quelccaya
  Ice Cap"; data updated 1992, NOAA template corrected 2015.
  URL:      https://www.ncei.noaa.gov/access/paleo-search/study/2551
  File:     https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/quelccaya/
            q83summ-noaa.txt
  Accessed: 2026-08-30
  Note:     This is reference (3) of the source figure, cited there to support
            the claim that the post-Maunder warming is "seen in many locations
            globally". It is included so that claim can be tested against the
            data rather than argued about.
            d18O at a tropical ice cap is NOT a clean temperature proxy -- the
            authors' own companion title is "A 1500-Year Record of Tropical
            Precipitation" -- so it is plotted on its own axis in permil, never
            on the temperature axis, and never converted to degrees.

recon_nh
recon_sh
  Hemispheric mean surface temperature reconstructions from real proxies,
  100-member ensembles, 1000-1999/2000 AD, anomalies wrt 1000-1999.
  Neukom et al. (2018) Nature Communications 9:5195, "Possible causes of data
  model discrepancy in the temperature history of the last Millennium".
  Proxy input is the PAGES 2k v2.0.0 database (PAGES2k Consortium 2017,
  Scientific Data 4:170088).
  URL:      https://www.ncei.noaa.gov/access/paleo-search/study/25455
  File:     https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/neukom2018/
            Real_proxy_recons/NH.txt  (and SH.txt)
  Accessed: 2026-08-30
  Note:     Semicolon separated, no header: column 1 is year CE, columns 2-101
            are ensemble members. The loader takes the per-year ensemble MEDIAN;
            the 5-95 percentile spread is computed separately and reported, so
            the uncertainty is shown rather than hidden behind a single line.
            These are hemispheric means. They are the comparison series for the
            claim that the post-Maunder excursion is global; HadCRUT5 begins in
            1850 and cannot address a 1680-1735 question.

sunspots
  SILSO version 2.0 yearly mean total sunspot number, 1700-present.
  Clette & Lefevre (2016) Solar Physics 291:2629-2651.
  URL:      https://www.sidc.be/SILSO/datafiles
  File:     https://www.sidc.be/SILSO/DATA/SN_y_tot_V2.0.txt
  Accessed: 2026-08-30
  Note:     Shown for context only. Not converted to a forcing and never drawn
            on a forcing or temperature axis.

==============================================================================
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA = Path(__file__).parent / "data"
SNAPSHOTS = DATA / "snapshots"

# Physical constants used to convert irradiance to forcing. These are geometry,
# not tuning: the Sun illuminates a cross-section pi*r^2 while the Earth
# radiates from 4*pi*r^2, and roughly 30% is reflected.
ALBEDO = 0.30
TSI_TO_FORCING = (1.0 - ALBEDO) / 4.0

# CO2 forcing coefficient, Myhre et al. (1998) simplified expression.
CO2_FORCING_COEFF = 5.35

# Single reference period for all forcing comparisons, per IPCC AR6 WG1 Ch.7.
FORCING_REF_YEAR = 1750

BASELINE = (1961, 1990)
XLIM = (1600, 2026)


class DataError(Exception):
    """Raised when a series fails validation. Never caught to keep a run alive."""


# ---------------------------------------------------------------------------
# Configuration: file names and the columns actually present in those files.
# Column names were read from the file headers, not guessed.
# ---------------------------------------------------------------------------

CFG = {
    "hadcet": {"file": "meantemp_monthly_totals.txt", "column": "Annual"},
    "hadcrut5": {"file": "hadcrut5_annual.csv", "column": "Anomaly (deg C)"},
    "uah_lt": {"file": "uah_lt_6.1.txt", "column": "Globe"},
    "co2_ice": {"file": "law2006.txt", "column": "CO2spl"},
    # co2_mlo, recon_* and sunspots are headerless; columns are positional and
    # the loader comment records which. None disables the header check.
    "co2_mlo": {"file": "co2_annmean_mlo.txt", "column": None},
    "tsi": {"file": "nrl2_tsi_P1Y.csv", "column": "irradiance (W/m^2)"},
    "quelccaya": {"file": "quelccaya_q83summ.txt", "column": "O18"},
    "recon_nh": {"file": "neukom2018_NH.txt", "column": None},
    "recon_sh": {"file": "neukom2018_SH.txt", "column": None},
    "sunspots": {"file": "SN_y_tot_V2.0.txt", "column": None},
}

# ---------------------------------------------------------------------------
# Expected coverage, content hash, physical range and minimum point count.
#
# Coverage alone is not enough: a wrong file spanning the right years passes a
# year check in silence. sha256 pins the exact bytes, value_range catches unit
# errors and column mix-ups, n_min catches a file that parsed into mostly NaN.
# ---------------------------------------------------------------------------

EXPECT = {
    "hadcet": {
        "first_year": 1659,
        "last_year_min": 2024,
        "sha256": "e0a74c142e73bf3d496246bce320710bb6781de301837da1fa0c3786326d471e",
        "value_range": (6.0, 13.0),  # absolute degC, English Midlands annual mean
        "n_min": 360,
    },
    "hadcrut5": {
        "first_year": 1850,
        "last_year_min": 2024,
        "sha256": "a1441ab5aef8f3baf43cca417bca271f6d674d4f235ea01129f0432c948882e0",
        "value_range": (-1.0, 2.0),
        "n_min": 170,
    },
    "uah_lt": {
        "first_year": 1979,
        "last_year_min": 2024,
        "sha256": "aaa1f4dd788807e706b720d3a00b783ec87604173c7ff73608f6bf49d23e0a51",
        "value_range": (-1.0, 1.5),
        "n_min": 44,
    },
    "co2_ice": {
        "first_year": 1600,
        "last_year_min": 2000,
        "sha256": "19ac7863487268e9ebf046cfc0c4fc481858469d82b6388ceee786181b1b5f05",
        "value_range": (270.0, 380.0),
        "n_min": 400,
    },
    "co2_mlo": {
        "first_year": 1959,
        "last_year_min": 2024,
        "sha256": "d84833898a9131921065912a9a33a8e1c4fb5791b5cd4f77ba5f71761cd717e5",
        "value_range": (310.0, 440.0),
        "n_min": 65,
    },
    "tsi": {
        "first_year": 1610,
        "last_year_min": 2020,
        "sha256": "48bada9d40d3a86d21b87f4356defde899fe73825bd64ebf0a34f945ed16e16d",
        "value_range": (1358.0, 1364.0),
        "n_min": 400,
    },
    "quelccaya": {
        "first_year": 1600,
        "last_year_min": 1980,
        "sha256": "feb80d6f979e2b9ec7cb4dfcb14c6353bf3dd2dae627d6ae8af1c0e53c8d67bc",
        "value_range": (-30.0, -5.0),  # permil VSMOW, tropical high-altitude ice
        "n_min": 350,
    },
    "recon_nh": {
        "first_year": 1600,
        "last_year_min": 1990,
        "sha256": "74bca8db0ae2cd9193f34f8da1f21c12ee2aa8745e4059d1216075dfbf082ea1",
        "value_range": (-1.5, 1.0),
        "n_min": 380,
    },
    "recon_sh": {
        "first_year": 1600,
        "last_year_min": 1990,
        "sha256": "25f6b96c0ea110c69b21b31c627a1ccd1894466a999b894a9acd8577ee36feab",
        "value_range": (-1.5, 1.0),
        "n_min": 380,
    },
    "sunspots": {
        "first_year": 1700,
        "last_year_min": 2024,
        "sha256": "7093eb8f5582a2d8e76580f32e6311ec758f0c1b969cadf99bcaba9f99726353",
        "value_range": (0.0, 300.0),
        "n_min": 300,
    },
}


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve(key: str) -> tuple[Path, str | None]:
    """Return the file to read and a staleness note if a snapshot was used.

    A normal run reads the working copy in data/. If it is missing, fall back to
    the newest pinned snapshot and report it. If neither exists, stop.
    """
    live = DATA / CFG[key]["file"]
    if live.exists():
        return live, None

    snap_dir = SNAPSHOTS / key
    if snap_dir.is_dir():
        dated = sorted(d for d in snap_dir.iterdir() if d.is_dir())
        if dated:
            cand = dated[-1] / CFG[key]["file"]
            if cand.exists():
                return cand, dated[-1].name

    raise DataError(
        f"{key}: no working copy at {live} and no snapshot under {snap_dir}. "
        f"Download it from the URL in the MANIFEST. Do not substitute another dataset."
    )


def check_column(key: str, path: Path) -> None:
    """Assert the column CFG names is really in the file's header.

    Providers rename and reorder columns between releases. Without this the CFG
    entry is documentation that can drift out of step with the loader, which is
    the failure mode that puts the wrong series in the figure under the right
    label.
    """
    col = CFG[key]["column"]
    if col is None:  # headerless file; the loader reads by position
        return
    with open(path, encoding="utf-8", errors="replace") as fh:
        head = "".join(next(fh, "") for _ in range(250))
    if col not in head:
        raise DataError(
            f"{key}: column {col!r} not found in the header of {path.name}. "
            f"The provider may have renamed or reordered columns. Read the header "
            f"and update CFG and the loader together -- do not guess."
        )


def check_manifest_completeness() -> None:
    """Every CFG series must be documented, and every documented series used."""
    doc = Path(__file__).read_text(encoding="utf-8")
    body = doc.split("MANIFEST", 1)[1] if "MANIFEST" in doc else ""
    missing = [k for k in CFG if f"\n{k}\n" not in body]
    if missing:
        raise DataError(
            f"MANIFEST has no entry for: {', '.join(sorted(missing))}. "
            f"A series without a citation and URL does not go in the figure."
        )
    unexpected = sorted(set(CFG) - set(EXPECT))
    if unexpected:
        raise DataError(f"CFG series with no EXPECT entry: {', '.join(unexpected)}")
    orphan = sorted(set(EXPECT) - set(CFG))
    if orphan:
        raise DataError(f"EXPECT entries with no CFG series: {', '.join(orphan)}")


# ---------------------------------------------------------------------------
# Loaders. Deliberately dumb: read file, select column, index by year.
# Cleverness in a loader hides data problems.
# ---------------------------------------------------------------------------


def load_hadcet(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=r"\s+", skiprows=4)
    s = pd.Series(df["Annual"].values, index=df["Year"].astype(int).values, dtype=float)
    return s[s > -99.0]  # -99.9 marks an incomplete year


def load_hadcrut5(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    return pd.Series(
        df["Anomaly (deg C)"].values, index=df["Time"].astype(int).values, dtype=float
    )


def load_uah_lt(path: Path) -> pd.Series:
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.split()
            if len(p) > 3 and p[0].isdigit() and len(p[0]) == 4 and p[1].isdigit():
                rows.append((int(p[0]), float(p[2])))
    if not rows:
        raise DataError("uah_lt: no monthly rows parsed; file format changed.")
    m = pd.DataFrame(rows, columns=["year", "anom"])
    counts = m.groupby("year")["anom"].count()
    annual = m.groupby("year")["anom"].mean()
    return annual[counts == 12]  # drop partial years rather than interpolate


def load_co2_ice(path: Path) -> pd.Series:
    rows = []
    started = False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("YearAD"):
                started = True
                continue
            if not started:
                continue
            p = line.split()
            if len(p) < 7:
                if rows:
                    break
                continue
            try:
                rows.append((int(float(p[4])), float(p[5])))
            except ValueError:
                break
    if not rows:
        raise DataError("co2_ice: spline table not found; file layout changed.")
    df = pd.DataFrame(rows, columns=["year", "co2"]).drop_duplicates("year")
    return pd.Series(df["co2"].values, index=df["year"].values, dtype=float)


def load_co2_mlo(path: Path) -> pd.Series:
    df = pd.read_csv(path, comment="#", sep=r"\s+", names=["year", "mean", "unc"])
    return pd.Series(df["mean"].values, index=df["year"].astype(int).values, dtype=float)


def load_tsi(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    days = df.iloc[:, 0].astype(float)
    # Header states "days since 1610-01-01"; values are mid-year points.
    year = (1610 + days / 365.2425).round().astype(int)
    s = pd.Series(df.iloc[:, 1].values, index=year.values, dtype=float)
    return s[~s.index.duplicated(keep="first")]


def load_quelccaya(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="	", comment="#")
    yr = pd.to_numeric(df["age_CE"], errors="coerce")
    d18 = pd.to_numeric(df["O18"], errors="coerce")
    ok = yr.notna() & d18.notna()
    s = pd.Series(d18[ok].values, index=yr[ok].astype(int).values, dtype=float)
    return s[~s.index.duplicated(keep="first")].sort_index()


def _load_recon_ensemble(path: Path) -> pd.DataFrame:
    """Year-indexed frame of the 100 ensemble members. Semicolon separated, no header."""
    d = pd.read_csv(path, sep=";", header=None)
    if d.shape[1] < 50:
        raise DataError(f"{path.name}: {d.shape[1]} columns, expected year + ~100 members.")
    return pd.DataFrame(
        d.iloc[:, 1:].astype(float).values, index=d.iloc[:, 0].astype(int).values
    )


def load_recon_nh(path: Path) -> pd.Series:
    return _load_recon_ensemble(path).median(axis=1)


def load_recon_sh(path: Path) -> pd.Series:
    return _load_recon_ensemble(path).median(axis=1)


def load_sunspots(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1], names=["y", "sn"])
    year = df["y"].astype(float).astype(int)
    s = pd.Series(df["sn"].values, index=year.values, dtype=float)
    return s[s >= 0]


LOADERS = {
    "hadcet": load_hadcet,
    "hadcrut5": load_hadcrut5,
    "uah_lt": load_uah_lt,
    "co2_ice": load_co2_ice,
    "co2_mlo": load_co2_mlo,
    "tsi": load_tsi,
    "quelccaya": load_quelccaya,
    "recon_nh": load_recon_nh,
    "recon_sh": load_recon_sh,
    "sunspots": load_sunspots,
}


# ---------------------------------------------------------------------------
# Validation gate
# ---------------------------------------------------------------------------


def validate(key: str, s: pd.Series, path: Path, digest: str) -> None:
    exp = EXPECT[key]

    if digest != exp["sha256"]:
        raise DataError(
            f"{key}: sha256 mismatch.\n"
            f"  expected {exp['sha256']}\n"
            f"  got      {digest}\n"
            f"  file     {path}\n"
            f"The provider may have reissued the dataset. Establish what changed "
            f"before updating the hash. Do not update it to silence this check."
        )

    s = s.dropna()
    if s.empty:
        raise DataError(f"{key}: no finite values after loading.")

    first, last = int(s.index.min()), int(s.index.max())
    if first > exp["first_year"]:
        raise DataError(
            f"{key}: starts {first}, expected coverage from {exp['first_year']}."
        )
    if last < exp["last_year_min"]:
        raise DataError(
            f"{key}: ends {last}, expected at least {exp['last_year_min']}."
        )

    lo, hi = exp["value_range"]
    out = s[(s < lo) | (s > hi)]
    if not out.empty:
        raise DataError(
            f"{key}: {len(out)} value(s) outside the physical range [{lo}, {hi}]; "
            f"first offender {out.index[0]}={out.iloc[0]:.4g}. "
            f"Likely a unit error or the wrong column."
        )

    if len(s) < exp["n_min"]:
        raise DataError(
            f"{key}: {len(s)} finite points, expected at least {exp['n_min']}."
        )


def load_all() -> tuple[dict[str, pd.Series], list[str]]:
    check_manifest_completeness()
    series: dict[str, pd.Series] = {}
    qc_rows: list[str] = []
    stale: list[str] = []

    for key in CFG:
        path, snap_date = resolve(key)
        digest = sha256_of(path)
        check_column(key, path)
        s = LOADERS[key](path).dropna().sort_index()
        s = s[(s.index >= XLIM[0] - 60) & (s.index <= XLIM[1])]
        validate(key, s, path, digest)
        series[key] = s
        if snap_date:
            stale.append(f"{key} (snapshot {snap_date})")
        qc_rows.append(
            "  {tag}{k:<10s} {fy:>5d}={fv:>10.4f}  {ly:>5d}={lv:>10.4f}  "
            "n={n:>5d}  sha={h}".format(
                tag="STALE: " if snap_date else "       ",
                k=key,
                fy=int(s.index.min()),
                fv=float(s.iloc[0]),
                ly=int(s.index.max()),
                lv=float(s.iloc[-1]),
                n=len(s),
                h=digest[:12],
            )
        )

    return series, qc_rows, stale


def main() -> int:
    try:
        series, qc_rows, stale = load_all()
    except DataError as exc:
        print(f"\nDataError: {exc}\n", file=sys.stderr)
        return 1

    print("\n" + "=" * 96)
    print("QC TABLE  (first year=value, last year=value, finite points, source sha256)")
    print("=" * 96)
    for row in qc_rows:
        print(row)
    print("=" * 96)
    if stale:
        print("STALE sources used: " + ", ".join(stale))
    print()

    try:
        d = derive(series)
    except DataError as exc:
        print(f"\nDataError: {exc}\n", file=sys.stderr)
        return 1

    print("DERIVED VALUES (every number annotated on the figure comes from here)")
    print("-" * 96)
    print(f"  CO2 forcing 1750->today            {d['f_co2_now']:+8.4f} W/m2")
    print(f"  Solar forcing 1750->today          {d['f_sol_now']:+8.4f} W/m2   "
          f"ratio {d['ratio']:.1f} : 1")
    print(f"  Solar forcing {d['trough_year']}->today (best case) "
          f"{d['f_sol_best']:+8.4f} W/m2   ratio {d['ratio_best']:.1f} : 1")
    print(f"  TSI trend since 1980               {d['tr_tsi'][0]:+8.4f} W/m2/decade "
          f"({d['tr_tsi'][1]}-{d['tr_tsi'][2]})")
    print(f"  HadCRUT5 trend since 1980          {d['tr_had'][0]:+8.4f} K/decade "
          f"({d['tr_had'][1]}-{d['tr_had'][2]})")
    print(f"  UAH trend since 1980               {d['tr_uah'][0]:+8.4f} K/decade "
          f"({d['tr_uah'][1]}-{d['tr_uah'][2]})")
    print(f"  Strongest 11-yr TSI                {d['tsi_peak']:.4f} in {d['tsi_peak_year']}; "
          f"{d['tsi_last_year']} ranks {d['tsi_rank']} of {d['tsi_n']}")
    print(f"  1680-1700 -> 1725-1735   HadCET    {d['exc']['cet']:+8.4f} K")
    print(f"                           NH recon  {d['exc']['nh']:+8.4f} K   "
          f"ratio {d['exc_ratio']:.1f} : 1")
    print(f"                           SH recon  {d['exc']['sh']:+8.4f} K")
    print(f"                           Quelccaya {d['exc_quel']:+8.4f} permil (NOT a temperature)")
    print("-" * 96 + "\n")

    out = Path(__file__).parent / "figure_climate_1600_2025.png"
    make_figure(d, out)
    print(f"Wrote {out}\n")
    return 0



# ---------------------------------------------------------------------------
# Derived quantities. Every number annotated on the figure comes from here,
# computed in the same run. Nothing below is typed in as a literal.
# ---------------------------------------------------------------------------


def cycle_average(s: pd.Series, window: int = 11) -> pd.Series:
    """Centred 11-year mean: one Schwabe cycle.

    The 11-year cycle is the largest amplitude in the TSI record and the
    smallest contributor to multidecadal temperature -- ocean thermal inertia
    damps it almost entirely. Comparisons against temperature or CO2 use this
    series, never the raw one.
    """
    return s.rolling(window, center=True, min_periods=window).mean().dropna()


def to_anomaly(s: pd.Series, label: str) -> pd.Series:
    """Re-reference to 1961-1990, warning when it happens."""
    lo, hi = BASELINE
    ref = s.loc[lo:hi]
    if len(ref) < (hi - lo + 1) * 0.9:
        raise DataError(
            f"{label}: only {len(ref)} years inside {lo}-{hi}; cannot re-reference. "
            f"Do not compare anomalies from different baselines."
        )
    print(f"  re-referencing {label} to {lo}-{hi} (offset {-ref.mean():+.4f})")
    return s - ref.mean()


def derive(series: dict[str, pd.Series]) -> dict:
    d: dict = {}

    tsi_raw = series["tsi"]
    tsi = cycle_average(tsi_raw)
    d["tsi_raw"], d["tsi"] = tsi_raw, tsi

    # --- forcings, single reference period, AR6 WG1 Ch.7 ---
    co2 = pd.concat([series["co2_ice"][series["co2_ice"].index < series["co2_mlo"].index.min()],
                     series["co2_mlo"]])
    d["co2"], d["co2_ice"], d["co2_mlo"] = co2, series["co2_ice"], series["co2_mlo"]

    c_ref = float(series["co2_ice"].loc[FORCING_REF_YEAR])
    t_ref = float(tsi.loc[FORCING_REF_YEAR])
    d["f_co2"] = CO2_FORCING_COEFF * np.log(co2 / c_ref)
    d["splice_year"] = int(series["co2_mlo"].index.min())
    d["f_co2_ice"] = CO2_FORCING_COEFF * np.log(series["co2_ice"] / c_ref)
    d["f_sol"] = (tsi - t_ref) * TSI_TO_FORCING

    d["f_co2_now"] = float(d["f_co2"].iloc[-1])
    d["f_sol_now"] = float(d["f_sol"].iloc[-1])
    d["ratio"] = d["f_co2_now"] / abs(d["f_sol_now"])

    # most solar-favourable framing: measure the Sun from its Maunder trough
    trough_y = int(tsi.loc[1600:1720].idxmin())
    d["trough_year"] = trough_y
    d["f_sol_best"] = float((tsi.iloc[-1] - tsi.loc[trough_y]) * TSI_TO_FORCING)
    d["ratio_best"] = d["f_co2_now"] / d["f_sol_best"]

    # --- temperatures on a common 1961-1990 baseline ---
    d["hadcrut5"] = series["hadcrut5"]
    d["cet"] = to_anomaly(series["hadcet"], "hadcet")
    d["nh"] = to_anomaly(series["recon_nh"], "recon_nh")
    d["sh"] = to_anomaly(series["recon_sh"], "recon_sh")

    # UAH begins 1979 and has no overlap with 1961-1990, so it cannot be
    # re-referenced. Offset it onto HadCRUT5 over their common years instead,
    # and report the offset rather than hiding it.
    uah = series["uah_lt"]
    common = uah.index.intersection(d["hadcrut5"].index)
    d["uah_offset"] = float(d["hadcrut5"].loc[common].mean() - uah.loc[common].mean())
    d["uah"] = uah + d["uah_offset"]
    d["uah_common"] = (int(common.min()), int(common.max()))

    d["quelccaya"] = series["quelccaya"]
    d["sunspots"] = series["sunspots"]

    # --- post-1980 divergence ---
    def slope(s: pd.Series, y0: int = 1980) -> tuple[float, int, int]:
        w = s.loc[y0:]
        b = float(np.polyfit(w.index.values.astype(float), w.values, 1)[0])
        return b * 10.0, int(w.index.min()), int(w.index.max())

    d["tr_tsi"] = slope(tsi)
    d["tr_fsol"] = slope(d["f_sol"])  # same trend in forcing units
    d["tr_had"] = slope(d["hadcrut5"])
    d["tr_uah"] = slope(d["uah"])

    # --- was the modern Sun the strongest on record? ---
    d["tsi_peak_year"] = int(tsi.idxmax())
    d["tsi_peak"] = float(tsi.max())
    d["tsi_last_year"] = int(tsi.index.max())
    d["tsi_last"] = float(tsi.iloc[-1])
    d["tsi_rank"] = int((tsi > tsi.iloc[-1]).sum()) + 1
    d["tsi_n"] = int(len(tsi))

    # --- the 1730 excursion, identical windows, per series ---
    def exc(s: pd.Series) -> float:
        return float(s.loc[1725:1735].mean() - s.loc[1680:1700].mean())

    d["exc"] = {"cet": exc(d["cet"]), "nh": exc(d["nh"]), "sh": exc(d["sh"])}
    d["exc_ratio"] = d["exc"]["cet"] / d["exc"]["nh"]
    d["exc_quel"] = exc(series["quelccaya"])  # permil, NOT a temperature

    return d


# ---------------------------------------------------------------------------
# Figure. Same architecture as the source figure -- one wide panel, 1600-2025,
# era shading, annotated turning points -- so the two can be read side by side.
# The corrections are in the axes, not the layout.
# ---------------------------------------------------------------------------

C_SOL, C_CO2 = "#E8890C", "#1B8A3A"
C_CET, C_GLB, C_REC, C_SAT = "#8B2FC9", "#111111", "#1F6FB8", "#C4243C"

# Real extent of the solar minima. Drawn from x-axis coordinates, never nudged.
MINIMA = [("Maunder Minimum", 1645, 1715), ("Dalton Min.", 1795, 1815)]

EXC_LABEL = {
    "cet": "HadCET\n(Central England)",
    "nh": "Neukom 2018\nN. Hemisphere",
    "sh": "Neukom 2018\nS. Hemisphere",
}


def make_figure(d: dict, out: Path) -> None:
    fig = plt.figure(figsize=(17.5, 11.6))
    gs = fig.add_gridspec(
        2, 4, height_ratios=[2.30, 1.0], hspace=0.60, wspace=0.34,
        left=0.055, right=0.945, top=0.850, bottom=0.075,
    )
    ax = fig.add_subplot(gs[0, :])
    axT = ax.twinx()

    exc_cet = d["exc"]["cet"]
    exc_nh = d["exc"]["nh"]

    # ---- era shading, drawn from x-axis coordinates ----
    for name, y0, y1 in MINIMA:
        ax.axvspan(y0, y1, color="#6A5ACD", alpha=0.10, zorder=0)
        ax.text((y0 + y1) / 2, 0.022, name + "\n" + str(y0) + "-" + str(y1),
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=8.5, color="#4A3F9F")

    # ---- forcings: left axis, both series in W/m2 on ONE scale ----
    lbl_co2 = "CO$_2$ forcing (ice core $\\to$ Mauna Loa)"
    lbl_sol = "Solar forcing, 11-yr mean"
    ax.plot(d["f_co2"].index, d["f_co2"].values, color=C_CO2, lw=2.6, label=lbl_co2)
    ov = d["f_co2_ice"].loc[d["splice_year"]:]
    ax.plot(ov.index, ov.values, color=C_CO2, lw=1.4, ls=(0, (5, 2)), alpha=0.85,
            label="Law Dome ice core to %d (overlap)" % int(ov.index.max()))
    ax.axvline(d["splice_year"], color=C_CO2, lw=0.9, ls=":", alpha=0.8)
    ax.annotate("splice: ice core to Mauna Loa, %d" % d["splice_year"],
                xy=(d["splice_year"], float(d["f_co2"].loc[d["splice_year"]])),
                xytext=(d["splice_year"] - 96, 0.92), fontsize=8.2, color=C_CO2,
                arrowprops=dict(arrowstyle="->", color=C_CO2, lw=1.0),
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=C_CO2, alpha=0.9))
    ax.plot(d["f_sol"].index, d["f_sol"].values, color=C_SOL, lw=2.6, label=lbl_sol)
    fr = (d["tsi_raw"] - float(d["tsi"].loc[FORCING_REF_YEAR])) * TSI_TO_FORCING
    # Line weight and opacity only -- the amplitude is untouched. The Schwabe
    # cycle is small in forcing terms and must stay that way; it just should not
    # be invisible.
    ax.plot(fr.index, fr.values, color="#C46A00", lw=1.15, alpha=0.75,
            label="Solar forcing, raw 11-yr cycle")
    ax.axhline(0, color="#999999", lw=0.7, ls=":", zorder=0)
    ax.set_ylabel("Radiative forcing (W/m$^2$)  —  both curves, one scale", fontsize=10.5)
    ax.set_ylim(-0.62, 3.10)
    ax.set_xlim(*XLIM)
    ax.set_xlabel("Year", fontsize=10.5)

    # ---- temperatures: right axis ----
    axT.plot(d["nh"].index, d["nh"].values, color=C_REC, lw=1.5, alpha=0.9,
             label="N. Hemisphere reconstruction (Neukom 2018)")
    axT.plot(d["cet"].index, d["cet"].values, color=C_CET, lw=0.9, alpha=0.75,
             label="HadCET - Central England ONLY, not the globe")
    axT.plot(d["hadcrut5"].index, d["hadcrut5"].values, color=C_GLB, lw=2.0,
             label="HadCRUT5 global mean surface temperature")
    axT.plot(d["uah"].index, d["uah"].values, color=C_SAT, lw=1.3, ls="--",
             label="UAH lower troposphere (offset %+.2f K)" % d["uah_offset"])
    axT.set_ylabel("Temperature anomaly (degC, wrt 1961-1990)", fontsize=10.5)
    axT.set_ylim(-2.7, 2.0)

    # ---- annotations. Every value AND every anchor comes from the loaded data.
    # An arrow tip is a claim about the year and value it sits over, so both
    # coordinates are read off the series, never typed in to look right.
    cet_peak_year = int(d["cet"].loc[1725:1735].idxmax())
    cet_peak_val = float(d["cet"].loc[cet_peak_year])
    axT.annotate(
        ("The 1730 peak is ENGLISH.\nHadCET %+.2f K   vs   N. Hemisphere %+.2f K\n"
         "= %.0fx larger in one region than hemispherically"
         % (exc_cet, exc_nh, d["exc_ratio"])),
        # xytext is in axT (temperature) coordinates, not ax (forcing) ones
        xy=(cet_peak_year, cet_peak_val), xytext=(1604, 1.93), fontsize=9.5, va="top",
        color=C_CET, arrowprops=dict(arrowstyle="->", color=C_CET, lw=1.3),
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C_CET, alpha=0.93))

    # Anchor on the solar-forcing curve at its own peak year, so the arrow sits
    # over the feature the text is about.
    pk_y = d["tsi_peak_year"]
    ax.annotate(
        ("Strongest Sun was %d, not today.\n%d ranks %dth of %d years.\n"
         "Since 1980: solar forcing %+.3f W/m$^2$/decade\n"
         "while HadCRUT5 %+.3f K/decade."
         % (pk_y, d["tsi_last_year"], d["tsi_rank"], d["tsi_n"],
            d["tr_fsol"][0], d["tr_had"][0])),
        xy=(pk_y, float(d["f_sol"].loc[pk_y])), xytext=(1776, 3.04),
        fontsize=9.5, va="top", color="#7A4A00",
        arrowprops=dict(arrowstyle="->", color="#7A4A00", lw=1.3,
                        connectionstyle="arc3,rad=-0.15"),
        bbox=dict(boxstyle="round,pad=0.4", fc="#FFF6E8", ec=C_SOL, alpha=0.95))

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axT.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", fontsize=8.6, framealpha=0.0,
              bbox_to_anchor=(0.5, -0.085), ncol=4, columnspacing=2.4,
              handlelength=2.6, borderpad=0.2)
    ax.text(0.5, -0.170,
            "All forcings referenced to %d (IPCC AR6 WG1 Ch.7). Temperatures are "
            "anomalies wrt 1961-1990. The raw solar cycle is drawn but not filled: "
            "the ocean damps it.\n"
            "UAH is lower troposphere, not a surface record; its record starts 1979, "
            "so it cannot be re-referenced to 1961-1990 and is offset onto HadCRUT5 "
            "over %d-%d."
            % (FORCING_REF_YEAR, d["uah_common"][0], d["uah_common"][1]),
            transform=ax.transAxes, ha="center", va="top", fontsize=7.6,
            color="#555555", linespacing=1.5)

    # ===================== bottom row =====================

    # (0) The Sun at its own scale.
    # The main panel shows solar and CO2 forcing on one axis, which is the
    # honest magnitude comparison but leaves the solar structure too small to
    # read. Rather than inflate that curve -- the source figure's error with
    # the sign reversed -- the same series gets its own axis here, with the
    # magnification stated so the two panels cannot be confused.
    a0 = fig.add_subplot(gs[1, 0])
    for _n, _y0, _y1 in MINIMA:
        a0.axvspan(_y0, _y1, color="#6A5ACD", alpha=0.12, zorder=0)
    a0.plot(fr.index, fr.values, color=C_SOL, lw=0.7, alpha=0.40)
    a0.plot(d["f_sol"].index, d["f_sol"].values, color=C_SOL, lw=2.2)
    a0.axhline(0, color="#999999", lw=0.7, ls=":", zorder=0)
    a0.set_xlim(*XLIM)
    a0.set_ylabel("Solar forcing (W/m$^2$)", color=C_SOL, fontsize=9)
    a0.tick_params(axis="y", colors=C_SOL, labelsize=8)
    a0.tick_params(axis="x", labelsize=8)
    _lo, _hi = a0.get_ylim()
    _mag = (ax.get_ylim()[1] - ax.get_ylim()[0]) / (_hi - _lo)
    a0.set_title("The Sun did vary - here it is at its own scale\n"
                 "same series as above, magnified %.0fx" % _mag,
                 fontsize=8.8, pad=6)
    a0.annotate("Maunder", xy=(1680, float(d["f_sol"].loc[1680])),
                xytext=(1618, _lo + 0.17 * (_hi - _lo)), fontsize=7.5,
                color="#4A3F9F",
                arrowprops=dict(arrowstyle="->", color="#4A3F9F", lw=0.9))
    a0.annotate("peak %d" % pk_y, xy=(pk_y, float(d["f_sol"].loc[pk_y])),
                xytext=(1793, _hi - 0.06 * (_hi - _lo)), fontsize=7.5,
                color="#7A4A00", ha="center",
                arrowprops=dict(arrowstyle="->", color="#7A4A00", lw=0.9))

    # (1) the divergence since 1980
    a1 = fig.add_subplot(gs[1, 1])
    a1t = a1.twinx()
    # Solar shown as FORCING here too. Raw TSI in W/m2 never shares a panel
    # with a climate quantity -- that conflation is one of the defects being
    # corrected, and it would be no better for being ours.
    ts, hd = d["f_sol"].loc[1980:], d["hadcrut5"].loc[1980:]
    a1.plot(ts.index, ts.values, color=C_SOL, lw=2.2)
    a1t.plot(hd.index, hd.values, color=C_GLB, lw=2.2)
    a1t.plot(d["uah"].loc[1980:].index, d["uah"].loc[1980:].values,
             color=C_SAT, lw=1.1, ls="--")
    for s_, ax_, c_ in ((ts, a1, C_SOL), (hd, a1t, C_GLB)):
        b_, a_ = np.polyfit(s_.index.values.astype(float), s_.values, 1)
        ax_.plot(s_.index, b_ * s_.index.values + a_, color=c_, lw=1.0, ls=":")
    a1.set_ylabel("Solar forcing, 11-yr mean (W/m$^2$)", color=C_SOL, fontsize=9)
    a1t.set_ylabel("Temp. anomaly (degC)", color=C_GLB, fontsize=9)
    a1.tick_params(axis="y", colors=C_SOL, labelsize=8)
    a1.tick_params(axis="x", labelsize=8)
    a1t.tick_params(axis="y", colors=C_GLB, labelsize=8)
    a1.set_title("Since 1980 the Sun declines while the Earth warms\n"
                 "solar %+.3f W/m$^2$/dec (to %d)\n"
                 "HadCRUT5 %+.3f  ·  UAH %+.3f K/dec"
                 % (d["tr_fsol"][0], d["tr_fsol"][2], d["tr_had"][0], d["tr_uah"][0]),
                 fontsize=8.8, pad=6)

    # (2) forcing magnitudes
    a2 = fig.add_subplot(gs[1, 2])
    names = ["CO$_2$\n(1750 to today)", "Solar\n(1750 to today)",
             "Solar, most\nfavourable case\n(%d to today)" % d["trough_year"]]
    vals = [d["f_co2_now"], d["f_sol_now"], d["f_sol_best"]]
    bars = a2.bar(names, vals, color=[C_CO2, C_SOL, C_SOL], width=0.62)
    bars[2].set_alpha(0.55)
    for b_, v in zip(bars, vals):
        a2.text(b_.get_x() + b_.get_width() / 2, v + 0.06, "%+.2f" % v,
                ha="center", fontsize=10, fontweight="bold")
    a2.set_ylabel("Radiative forcing (W/m$^2$)", fontsize=9)
    a2.set_ylim(0, max(vals) * 1.32)
    a2.tick_params(labelsize=8)
    a2.set_title("CO$_2$ forcing is %.0fx the solar forcing\n"
                 "(%.0fx even from the Maunder trough)"
                 % (d["ratio"], d["ratio_best"]), fontsize=9.5, pad=7)

    # (3) the 1730 excursion, identical windows, per series
    a3 = fig.add_subplot(gs[1, 3])
    keys = ["cet", "nh", "sh"]
    vs = [d["exc"][k] for k in keys]
    bars = a3.bar([EXC_LABEL[k] for k in keys], vs, color=[C_CET, C_REC, C_REC], width=0.6)
    for b_, v in zip(bars, vs):
        a3.text(b_.get_x() + b_.get_width() / 2, v + 0.035, "%+.3f K" % v,
                ha="center", fontsize=9.5, fontweight="bold")
    a3.set_ylabel("Warming, 1680-1700 to 1725-1735 (K)", fontsize=9)
    a3.set_ylim(0, max(vs) * 1.42)
    a3.tick_params(labelsize=8)
    a3.set_title("The same event, measured the same way\n"
                 "England warms %.0fx more than the hemisphere" % d["exc_ratio"],
                 fontsize=8.8, pad=6)
    a3.text(0.5, 0.52,
            "Quelccaya ice cap, Peru - the source figure's own\n"
            "reference (3) for \"seen in many locations globally\":\n"
            "d18O moves %+.2f permil over the same window,\n"
            "i.e. the opposite direction. It is also a\n"
            "precipitation proxy, not a thermometer."
            % d["exc_quel"],
            transform=a3.transAxes, ha="center", va="center", fontsize=6.8,
            bbox=dict(boxstyle="round,pad=0.4", fc="#FFF9E6", ec="#D9B24C"))

    for a_ in (ax, a0, a1, a2, a3):
        a_.grid(alpha=0.16, lw=0.6)
        a_.set_axisbelow(True)

    fig.text(0.055, 0.972, "Was it the Sun?", fontsize=21, fontweight="bold", va="top")
    fig.text(0.055, 0.934,
             "The same 1600-2025 chart as the version circulating on X, rebuilt with every "
             "series loaded from its published archive and every axis in comparable units.",
             fontsize=10.5, va="top", color="#333333")
    fig.text(0.055, 0.899,
             "Corrections: (1) raw TSI in W/m2 is not a radiative forcing - the conversion is "
             "dTSI/4x(1-albedo) ~ 0.175;   (2) Central England is not the globe;   "
             "(3) the Spoerer Minimum lies outside a 1600-start axis.",
             fontsize=9, va="top", color="#555555")
    fig.text(0.945, 0.972,
             "Sources: HadCET · HadCRUT5 · UAH v6.1 · Law Dome (MacFarling Meure 2006) · "
             "NOAA Mauna Loa\nNRLTSI2 (Coddington 2016) · SILSO · Quelccaya (Thompson 1986) · "
             "Neukom 2018 / PAGES 2k v2\nFull citations, URLs and SHA-256 of every file: "
             "MANIFEST in climate_figure.py",
             fontsize=7.8, va="top", ha="right", color="#444444", linespacing=1.6)

    fig.savefig(out, dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
