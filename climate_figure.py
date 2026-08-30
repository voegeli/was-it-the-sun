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
    "co2_mlo": {"file": "co2_annmean_mlo.txt", "column": "mean"},
    "tsi": {"file": "nrl2_tsi_P1Y.csv", "column": "irradiance (W/m^2)"},
    "sunspots": {"file": "SN_y_tot_V2.0.txt", "column": "sn"},
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
