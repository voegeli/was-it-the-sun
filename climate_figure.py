"""
Solar activity vs CO2 as drivers of observed warming, 1600-2025.

A corrected rebuild of a widely circulated figure ("Our Changing Climate", slide
16, posted to X by @bootcanyon on 2026-08-26) which argues that the warming since
the Little Ice Age is solar in origin. The original is not redistributed here;
its defects are quoted as text and tested against the data. This script keeps
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

The MANIFEST is the structured dict below, not prose: it is the single source
of citations, URLs and access dates, it is checked at startup, and fetch_data.py
imports its URLs so the downloader cannot drift away from the documentation.
"""

from __future__ import annotations

import hashlib
import sys
from datetime import date
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
# MANIFEST -- provenance, one entry per CFG series.
#
# Prose is not checkable. Every field here is asserted non-empty at startup by
# check_manifest_completeness(), and file_url is the ONE place a download URL is
# written: fetch_data.py imports these rather than keeping its own copy. That
# removes the class of bug where the documented URL and the fetched URL differ.
# ---------------------------------------------------------------------------

MANIFEST = {
    "hadcet": {
        "citation": "Manley (1974) Q.J.R.Meteorol.Soc. 100:389-405; Parker et al. "
                    "(1992) Int.J.Climatol. 12:317-342; Parker & Horton (2005) "
                    "Int.J.Climatol. 25:1173-1188.",
        "landing": "https://www.metoffice.gov.uk/hadobs/hadcet/",
        "file_url": "https://www.metoffice.gov.uk/hadobs/hadcet/data/"
                    "meantemp_monthly_totals.txt",
        "accessed": "2026-08-30",
        "sha256": "e0a74c142e73bf3d496246bce320710bb6781de301837da1fa0c3786326d471e",
        "note": "Mean Central England Temperature, absolute degC. REGIONAL "
                "(English Midlands), not a global record. Converted here to an "
                "anomaly on its own 1961-1990 mean.",
    },
    "hadcrut5": {
        "citation": "Morice et al. (2021) JGR Atmospheres 126:e2019JD032361.",
        "landing": "https://www.metoffice.gov.uk/hadobs/hadcrut5/",
        "file_url": "https://www.metoffice.gov.uk/hadobs/hadcrut5/data/"
                    "HadCRUT.5.0.2.0/analysis/diagnostics/"
                    "HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv",
        "accessed": "2026-08-30",
        "sha256": "a1441ab5aef8f3baf43cca417bca271f6d674d4f235ea01129f0432c948882e0",
        "note": "Global mean surface temperature anomaly, already on the "
                "1961-1990 baseline; no re-referencing applied. The source figure "
                "labels its global series GloSAT; GloSAT needs a free CEDA account "
                "and is not used here. This series is labelled HadCRUT5 because "
                "that is what it is.",
    },
    "uah_lt": {
        "citation": "Spencer, Christy & Braswell (2017) Asia-Pac.J.Atmos.Sci. "
                    "53:121-130.",
        "landing": "https://www.nsstc.uah.edu/climate/",
        "file_url": "https://www.nsstc.uah.edu/data/msu/v6.1/tlt/uahncdc_lt_6.1.txt",
        "accessed": "2026-08-30",
        "sha256": "aaa1f4dd788807e706b720d3a00b783ec87604173c7ff73608f6bf49d23e0a51",
        "note": "Satellite lower-troposphere anomaly v6.1, published on a "
                "1991-2020 baseline. It CANNOT be re-referenced to 1961-1990 -- the "
                "record starts 1978, so there is no overlap with the baseline. It "
                "is offset onto HadCRUT5 over their common years and that offset is "
                "printed in the QC table. Lower troposphere is not a surface "
                "record; drawn dashed and labelled as a different quantity.",
    },
    "co2_ice": {
        "citation": "MacFarling Meure et al. (2006) Geophys.Res.Lett. 33:L14810; "
                    "Etheridge et al. (1996) JGR 101:4115-4128.",
        "landing": "https://www.ncei.noaa.gov/access/paleo-search/study/9959",
        "file_url": "https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/"
                    "law/law2006.txt",
        "accessed": "2026-08-30",
        "sha256": "19ac7863487268e9ebf046cfc0c4fc481858469d82b6388ceee786181b1b5f05",
        "note": "Law Dome ice core CO2, published smoothing spline (Enting 1987 "
                "method). The spline is the authors' product, not a smoothing "
                "applied here. Used before the Mauna Loa record begins; the splice "
                "is drawn with both series visible in overlap.",
    },
    "co2_mlo": {
        "citation": "Keeling et al. (1976); NOAA GML, Lan, Tans & Thoning.",
        "landing": "https://gml.noaa.gov/ccgg/trends/",
        "file_url": "https://gml.noaa.gov/webdata/ccgg/trends/co2/"
                    "co2_annmean_mlo.txt",
        "accessed": "2026-08-30",
        "sha256": "d84833898a9131921065912a9a33a8e1c4fb5791b5cd4f77ba5f71761cd717e5",
        "note": "NOAA GML Mauna Loa annual mean CO2.",
    },
    "tsi": {
        "citation": "Coddington et al. (2016) Bull.Amer.Meteor.Soc. 97:1265-1282; "
                    "Lean (2000) Geophys.Res.Lett. 27:2425-2428.",
        "landing": "https://lasp.colorado.edu/lisird/data/nrl2_tsi_P1Y",
        "file_url": "https://lasp.colorado.edu/lisird/latis/dap/nrl2_tsi_P1Y.csv",
        "accessed": "2026-08-30",
        "sha256": "48bada9d40d3a86d21b87f4356defde899fe73825bd64ebf0a34f945ed16e16d",
        "note": "NRLTSI2 total solar irradiance reconstruction, annual. Time axis "
                "is days since 1610-01-01, converted to calendar year in the "
                "loader. Absolute scale ~1361 W/m2. This is a low-variability "
                "reconstruction, and that choice is contestable: higher-amplitude "
                "reconstructions of the Maunder-to-present change exist. "
                "robustness.py block 8 computes what the objection is worth -- the "
                "solar term would need a dTSI about 19x larger, near 1% of the "
                "solar constant, to match the CO2 forcing. The post-1980 sign "
                "result is LESS dependent on that choice than the amplitude "
                "comparison is -- NRLTSI2 is calibrated against direct space-based "
                "observations over the satellite era and reproduces their "
                "variability -- but it is not independent of it. NRLTSI2 remains a "
                "proxy model built on sunspot and facular indices regressed onto "
                "those observations, and this pipeline loads no alternative TSI "
                "composite, so the defensible claim is that the trend since 1980 is "
                "negative WITHIN NRLTSI2, not that it is model-independent. "
                "See https://lasp.colorado.edu/lisird/data/nrl2_files/ and "
                "Coddington et al. (2016).",
    },
    "quelccaya": {
        "citation": "Thompson, Mosley-Thompson, Dansgaard & Grootes (1986) Science "
                    "234:361-364; data updated 1992, NOAA template corrected 2015.",
        "landing": "https://www.ncei.noaa.gov/access/paleo-search/study/2551",
        "file_url": "https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/"
                    "quelccaya/q83summ-noaa.txt",
        "accessed": "2026-08-30",
        "sha256": "feb80d6f979e2b9ec7cb4dfcb14c6353bf3dd2dae627d6ae8af1c0e53c8d67bc",
        "note": "Quelccaya Ice Cap summit core, annual d18O, Peru 13.9S 5670 m. "
                "This is reference (3) of the source figure, cited there for the "
                "claim that the post-Maunder warming is 'seen in many locations "
                "globally', and is included so that claim can be tested rather than "
                "argued about. d18O at a tropical ice cap is NOT a temperature "
                "proxy -- the authors' companion title is 'A 1500-Year Record of "
                "Tropical Precipitation' -- so it is reported in permil on its own "
                "terms, never on the temperature axis, never converted to degrees, "
                "and its sign is NOT read as warming or cooling.",
    },
    "recon_nh": {
        "citation": "Neukom et al. (2018) Nature Communications 9:5195; proxy input "
                    "is PAGES 2k v2.0.0 (PAGES2k Consortium 2017, Sci.Data 4:170088).",
        "landing": "https://www.ncei.noaa.gov/access/paleo-search/study/25455",
        "file_url": "https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/"
                    "neukom2018/Real_proxy_recons/NH.txt",
        "accessed": "2026-08-30",
        "sha256": "74bca8db0ae2cd9193f34f8da1f21c12ee2aa8745e4059d1216075dfbf082ea1",
        "note": "Northern Hemisphere mean temperature, 100-member ensemble, "
                "anomalies wrt 1000-1999. Semicolon separated, no header: column 1 "
                "is year CE, columns 2-101 are members. The ensemble is kept, not "
                "collapsed: the median is plotted and the 5-95 percentile range of "
                "the 1680-1735 excursion is computed, printed in QC and drawn as an "
                "interval. HadCRUT5 begins in 1850 and cannot address a 1680-1735 "
                "question; this series can.",
    },
    "recon_sh": {
        "citation": "Neukom et al. (2018) Nature Communications 9:5195; proxy input "
                    "is PAGES 2k v2.0.0 (PAGES2k Consortium 2017, Sci.Data 4:170088).",
        "landing": "https://www.ncei.noaa.gov/access/paleo-search/study/25455",
        "file_url": "https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/"
                    "neukom2018/Real_proxy_recons/SH.txt",
        "accessed": "2026-08-30",
        "sha256": "25f6b96c0ea110c69b21b31c627a1ccd1894466a999b894a9acd8577ee36feab",
        "note": "Southern Hemisphere counterpart of recon_nh, same format and "
                "treatment. Its 1680-1735 excursion interval spans zero, so the "
                "sign of that excursion is NOT established and must not be reported "
                "as if it were.",
    },
}

# Deliberately NOT part of this figure:
#
#   SILSO sunspot number (Clette & Lefevre 2016, https://www.sidc.be/SILSO/).
#   It was loaded and validated in an earlier revision but never plotted or used
#   in any computation, while the MANIFEST claimed it was "shown for context".
#   Sunspot number is a proxy for the same solar activity NRLTSI2 already
#   reconstructs, it cannot be converted to a radiative forcing, and a second
#   solar curve on this figure would add visual weight without adding evidence.
#   It is removed rather than left half-registered -- a series that looks checked
#   but is unused is worse than an absent one. Restoring it means restoring its
#   CFG, MANIFEST, EXPECT and QC entries together, in a clearly separated context
#   panel, never on a forcing or temperature axis.


# ---------------------------------------------------------------------------
# Configuration: file names and the columns actually present in those files.
# Column names were read from the file headers, not guessed.
# ---------------------------------------------------------------------------

CFG = {
    "hadcet": {"file": "meantemp_monthly_totals.txt", "column": "Annual"},
    "hadcrut5": {"file": "hadcrut5_annual.csv", "column": "Anomaly (deg C)"},
    "uah_lt": {"file": "uah_lt_6.1.txt", "column": "Globe"},
    "co2_ice": {"file": "law2006.txt", "column": "CO2spl"},
    # co2_mlo and recon_* are headerless; columns are positional and
    # the loader comment records which. None disables the header check.
    "co2_mlo": {"file": "co2_annmean_mlo.txt", "column": None},
    "tsi": {"file": "nrl2_tsi_P1Y.csv", "column": "irradiance (W/m^2)"},
    "quelccaya": {"file": "quelccaya_q83summ.txt", "column": "O18"},
    "recon_nh": {"file": "neukom2018_NH.txt", "column": None},
    "recon_sh": {"file": "neukom2018_SH.txt", "column": None},
}

# ---------------------------------------------------------------------------
# Expected coverage, content hash, physical range and minimum point count.
#
# Coverage alone is not enough: a wrong file spanning the right years passes a
# year check in silence. The sha256 that pins the exact bytes lives in MANIFEST
# and is NOT copied here -- two hand-maintained copies of one hash drift, and a
# drifted hash is worse than none because it still looks checked. value_range
# catches unit errors and column mix-ups, n_min catches a file that parsed into
# mostly NaN.
# ---------------------------------------------------------------------------

EXPECT = {
    "hadcet": {
        "first_year": 1659,
        "last_year_min": 2024,
        "value_range": (6.0, 13.0),  # absolute degC, English Midlands annual mean
        "n_min": 360,
    },
    "hadcrut5": {
        "first_year": 1850,
        "last_year_min": 2024,
        "value_range": (-1.0, 2.0),
        "n_min": 170,
    },
    "uah_lt": {
        "first_year": 1979,
        "last_year_min": 2024,
        "value_range": (-1.0, 1.5),
        "n_min": 44,
    },
    "co2_ice": {
        "first_year": 1600,
        "last_year_min": 2000,
        "value_range": (270.0, 380.0),
        "n_min": 400,
    },
    "co2_mlo": {
        "first_year": 1959,
        "last_year_min": 2024,
        "value_range": (310.0, 440.0),
        "n_min": 65,
    },
    "tsi": {
        "first_year": 1610,
        "last_year_min": 2020,
        "value_range": (1358.0, 1364.0),
        "n_min": 400,
    },
    "quelccaya": {
        "first_year": 1600,
        "last_year_min": 1980,
        "value_range": (-30.0, -5.0),  # permil VSMOW, tropical high-altitude ice
        "n_min": 350,
    },
    "recon_nh": {
        "first_year": 1600,
        "last_year_min": 1990,
        "value_range": (-1.5, 1.0),
        "n_min": 380,
    },
    "recon_sh": {
        "first_year": 1600,
        "last_year_min": 1990,
        "value_range": (-1.5, 1.0),
        "n_min": 380,
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


REQUIRED_MANIFEST_FIELDS = (
    "citation", "landing", "file_url", "accessed", "sha256", "note",
)


def check_manifest_completeness() -> None:
    """Assert the MANIFEST really documents every series, field by field.

    The previous version only looked for the series name somewhere in the source
    text, which a comment mentioning the name would satisfy. That is not
    provenance. These checks fail loudly instead:

      - CFG, MANIFEST and EXPECT cover exactly the same set of series
      - every required field is present, a string, and not blank
      - URLs are URLs, and the file_url is not merely the landing page
      - the access date parses as a real ISO date
      - the sha256 is 64 lowercase hex characters, and nothing else
    """
    cfg, man, exp = set(CFG), set(MANIFEST), set(EXPECT)

    if cfg - man:
        raise DataError(
            f"MANIFEST has no entry for: {', '.join(sorted(cfg - man))}. "
            f"A series without a citation and URL does not go in the figure."
        )
    if man - cfg:
        raise DataError(
            f"MANIFEST documents series that CFG does not load: "
            f"{', '.join(sorted(man - cfg))}. Remove the orphan or restore the series."
        )
    if cfg - exp:
        raise DataError(f"CFG series with no EXPECT entry: {', '.join(sorted(cfg - exp))}")
    if exp - cfg:
        raise DataError(f"EXPECT entries with no CFG series: {', '.join(sorted(exp - cfg))}")

    for key in sorted(cfg):
        entry = MANIFEST[key]
        for field in REQUIRED_MANIFEST_FIELDS:
            if field not in entry:
                raise DataError(f"MANIFEST[{key!r}] is missing {field!r}.")
            val = entry[field]
            if not isinstance(val, str) or not val.strip():
                raise DataError(
                    f"MANIFEST[{key!r}][{field!r}] is empty. An empty citation or "
                    f"URL is not provenance; fill it in or drop the series."
                )

        for field in ("landing", "file_url"):
            if not entry[field].startswith(("http://", "https://")):
                raise DataError(
                    f"MANIFEST[{key!r}][{field!r}] is not a URL: {entry[field]!r}"
                )
        if entry["file_url"].rstrip("/") == entry["landing"].rstrip("/"):
            raise DataError(
                f"MANIFEST[{key!r}]: file_url is the landing page. fetch_data.py "
                f"downloads file_url, so it must point at the actual file."
            )

        try:
            date.fromisoformat(entry["accessed"])
        except ValueError as exc:
            raise DataError(
                f"MANIFEST[{key!r}]['accessed'] is not an ISO date: "
                f"{entry['accessed']!r}"
            ) from exc

        digest = entry["sha256"]
        if len(digest) != 64 or not all(c in "0123456789abcdef" for c in digest):
            raise DataError(
                f"MANIFEST[{key!r}]['sha256'] is not 64 lowercase hex characters: "
                f"{digest!r}. Provenance without a usable content hash only proves "
                f"a file was named correctly, not that it is the right file."
            )


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
}


# ---------------------------------------------------------------------------
# Validation gate
# ---------------------------------------------------------------------------


def validate(key: str, s: pd.Series, path: Path, digest: str) -> None:
    exp = EXPECT[key]

    want = MANIFEST[key]["sha256"]  # single source; EXPECT does not duplicate it
    if digest != want:
        raise DataError(
            f"{key}: sha256 mismatch.\n"
            f"  expected {want}\n"
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
    cy = d["compare_year"]
    print(f"  comparison end year (both terms)     {cy}   "
          f"[CO2 data to {int(d['f_co2'].index.max())}, "
          f"centred 11-yr solar mean to {int(d['f_sol'].index.max())}]")
    print(f"  CO2 forcing {FORCING_REF_YEAR}->{cy}              "
          f"{d['f_co2_now']:+8.4f} W/m2")
    print(f"  Solar forcing {FORCING_REF_YEAR}->{cy}            "
          f"{d['f_sol_now']:+8.4f} W/m2   ratio {d['ratio']:.1f} : 1")
    print(f"  Solar forcing {d['trough_year']}->{cy} (best case) "
          f"{d['f_sol_best']:+8.4f} W/m2   ratio {d['ratio_best']:.1f} : 1")
    print(f"  TSI trend since 1980               {d['tr_tsi'][0]:+8.4f} W/m2/decade "
          f"({d['tr_tsi'][1]}-{d['tr_tsi'][2]})")
    print(f"  HadCRUT5 trend since 1980          {d['tr_had'][0]:+8.4f} K/decade "
          f"({d['tr_had'][1]}-{d['tr_had'][2]})")
    print(f"  UAH trend since 1980               {d['tr_uah'][0]:+8.4f} K/decade "
          f"({d['tr_uah'][1]}-{d['tr_uah'][2]})")
    print(f"  Strongest 11-yr TSI                {d['tsi_peak']:.4f} in "
          f"{d['tsi_peak_year']}; {d['tsi_last_year']} ranks {d['tsi_rank']} "
          f"of {d['tsi_n']}")
    print(f"  1680-1700 -> 1725-1735   HadCET    {d['exc']['cet']:+8.4f} K "
          f"(single series, no ensemble)")
    for key, name in (("nh", "NH recon "), ("sh", "SH recon ")):
        ci = d["exc_ci"][key]
        print(f"                           {name} {ci['median']:+8.4f} K   "
              f"5-95% [{ci['p5']:+.4f}, {ci['p95']:+.4f}]  "
              f"{ci['frac_pos'] * 100:.0f}% of {ci['n']} members positive")
    print(f"                           across {d['alt_n']} alternative window "
          f"definitions: HadCET median {d['alt_cet'][0]:+.3f} K "
          f"[{d['alt_cet'][1]:+.3f}, {d['alt_cet'][2]:+.3f}], "
          f"NH median {d['alt_nh'][0]:+.3f} K "
          f"[{d['alt_nh'][1]:+.3f}, {d['alt_nh'][2]:+.3f}]")
    print("                           the published window is at the favourable "
          "end for HadCET; the NH result is not sensitive to it")
    print(f"                           HadCET excursion is the {d['cet_pct']:.0f}th "
          f"percentile of {d['cet_windows_n']} CET 55-year windows")
    print(f"                           ratio     HadCET / NH median "
          f"{d['exc_ratio']:.1f} : 1;  vs NH 95th pct "
          f"{d['exc_ratio_min']:.1f} : 1")
    print("                           NOTE: the ratio is reported but NOT used as "
          "the headline -- its\n                           denominator is near zero "
          "in many windows, so the quotient is unstable.")
    if d["exc_ci"]["sh"]["p5"] < 0 < d["exc_ci"]["sh"]["p95"]:
        print("                           NOTE: the SH interval spans zero -- the "
              "SIGN of the SH excursion is not established.")
    print(f"                           Quelccaya {d['exc_quel']:+8.4f} permil "
          f"(d18O; not a temperature, sign not interpreted)")
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


# ---------------------------------------------------------------------------
# The "55-year excursion window", defined once.
#
# For a start year y0:
#     baseline = mean of the series over  y0      .. y0 + 20   (21 years)
#     peak     = mean of the series over  y0 + 45 .. y0 + 55   (11 years)
#     excursion = peak - baseline
#
# The published comparison is the y0 = 1680 member of this family: 1680-1700
# against 1725-1735. Everything that quotes a percentile, a distribution or a
# window count must go through these functions, so the figure and audit.py
# cannot disagree about which windows exist.
# ---------------------------------------------------------------------------

EXC_BASELINE = (0, 20)
EXC_PEAK = (45, 55)
EXC_LEN = 55
EXC_START = 1680  # the window the source figure's claim is about


def excursion(s: pd.Series, y0: int = EXC_START) -> float:
    """Excursion of one series over the window starting at y0."""
    base = s.loc[y0 + EXC_BASELINE[0]:y0 + EXC_BASELINE[1]].mean()
    peak = s.loc[y0 + EXC_PEAK[0]:y0 + EXC_PEAK[1]].mean()
    return float(peak - base)


def excursion_windows(s: pd.Series) -> pd.Series:
    """Every window the series supports, indexed by start year.

    The start year runs from the first year of the record; there is no
    hand-picked lower bound, because a hand-picked bound is how the figure and
    an audit of the figure end up counting different populations.
    """
    out = {}
    for y0 in range(int(s.index.min()), int(s.index.max()) - EXC_LEN):
        base = s.loc[y0 + EXC_BASELINE[0]:y0 + EXC_BASELINE[1]].mean()
        peak = s.loc[y0 + EXC_PEAK[0]:y0 + EXC_PEAK[1]].mean()
        if np.isfinite(base) and np.isfinite(peak):
            out[y0] = float(peak - base)
    return pd.Series(out, dtype=float)


def derive(series: dict[str, pd.Series]) -> dict:
    d: dict = {"exc_ci": {}}

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

    # The two series do not end in the same year: CO2 runs to the present while
    # the CENTRED 11-year solar mean necessarily stops half a cycle short of the
    # end of the TSI record. Comparing a 2025 CO2 value with a 2018 solar value
    # and calling both "today" overstates the ratio. Both terms are therefore
    # evaluated at the last year BOTH are defined, and that year is what the
    # figure prints -- never the word "today".
    #
    # The alternative, extending the solar mean with a trailing or partial
    # window, is deliberately not taken: it would change what the solar number
    # means without changing its label.
    d["compare_year"] = int(min(d["f_co2"].index.max(), d["f_sol"].index.max()))
    cy = d["compare_year"]
    if cy < d["f_co2"].index.max():
        print(f"  forcing comparison pinned to {cy}: CO2 runs to "
              f"{int(d['f_co2'].index.max())} but the centred 11-yr solar mean "
              f"ends {int(d['f_sol'].index.max())}")

    d["f_co2_now"] = float(d["f_co2"].loc[cy])
    d["f_sol_now"] = float(d["f_sol"].loc[cy])
    d["ratio"] = d["f_co2_now"] / abs(d["f_sol_now"])

    # most solar-favourable framing: measure the Sun from its Maunder trough,
    # still at the common end year
    trough_y = int(tsi.loc[1600:1720].idxmin())
    d["trough_year"] = trough_y
    d["f_sol_best"] = float((tsi.loc[cy] - tsi.loc[trough_y]) * TSI_TO_FORCING)
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
    exc = excursion  # the y0 = EXC_START member of the shared window family

    d["exc"] = {"cet": exc(d["cet"]), "nh": exc(d["nh"]), "sh": exc(d["sh"])}

    # The reconstructions are 100-member ensembles. Collapsing them to a median
    # and quoting one number hides the only thing that says how well the
    # excursion is resolved, so the spread is computed from the same files in
    # the same run and carried through to the QC table and the figure.
    for key in ("nh", "sh"):
        path, _ = resolve("recon_" + key)
        ens = _load_recon_ensemble(path)
        a = ens.loc[1680:1700].mean(axis=0)
        b = ens.loc[1725:1735].mean(axis=0)
        diff = (b - a).to_numpy(dtype=float)
        d["exc_ci"][key] = {
            "median": float(np.median(diff)),
            "p5": float(np.percentile(diff, 5)),
            "p95": float(np.percentile(diff, 95)),
            "frac_pos": float((diff > 0).mean()),
            "n": int(diff.size),
        }

    # The bar and its interval must be the same statistic. d["exc"]["nh"] is the
    # excursion OF the median series; exc_ci["nh"]["median"] is the median OF the
    # per-member excursions. They differ slightly, so the ensemble statistic is
    # used for both the plotted bar and the interval, and for the ratio.
    for key in ("nh", "sh"):
        d["exc"][key] = d["exc_ci"][key]["median"]

    # The RATIO of the two excursions is deliberately not the headline. Its
    # denominator is a hemispheric mean that is near zero in many windows, so the
    # quotient is unstable: across sliding 55-year windows its distribution runs
    # from below 1 to above 25, and a large value is not by itself evidence that
    # 1730 was special. It is computed and reported, but the figure leads with
    # the two absolute numbers instead.
    d["exc_ratio"] = d["exc"]["cet"] / d["exc"]["nh"]
    d["exc_ratio_min"] = d["exc"]["cet"] / d["exc_ci"]["nh"]["p95"]

    # The stable statement is about the SIZE of the English excursion: where does
    # it sit among every 55-year window the CET record offers, measured exactly
    # the same way? That needs no division and no second series.
    # How much of the published excursion is the window choice? Vary the start
    # year, the baseline length and the peak offsets, and carry the spread into
    # the figure. The published pair sits at the favourable end for England, and
    # saying so is cheaper than being caught not saying it.
    alt_cet, alt_nh = [], []
    for y0 in range(1670, 1691, 5):
        for base_len in (15, 20, 25):
            for lo, hi in ((40, 50), (45, 55), (50, 60)):
                for key, out in (("cet", alt_cet), ("nh", alt_nh)):
                    s = d[key]
                    out.append(float(s.loc[y0 + lo:y0 + hi].mean()
                                     - s.loc[y0:y0 + base_len].mean()))
    d["alt_n"] = len(alt_cet)
    d["alt_cet"] = (float(np.median(alt_cet)), float(min(alt_cet)), float(max(alt_cet)))
    d["alt_nh"] = (float(np.median(alt_nh)), float(min(alt_nh)), float(max(alt_nh)))

    win = excursion_windows(d["cet"])
    d["cet_windows"] = win
    d["cet_windows_n"] = int(win.size)
    d["cet_windows_first"] = int(win.index.min())
    d["cet_windows_last"] = int(win.index.max())
    d["cet_pct"] = float((win.to_numpy() < d["exc"]["cet"]).mean() * 100.0)
    # d18O in permil. Its sign is NOT converted to a temperature direction: no
    # proxy model for tropical d18O is loaded here, and the authors describe the
    # record as precipitation-influenced.
    d["exc_quel"] = exc(series["quelccaya"])

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
        left=0.055, right=0.945, top=0.850, bottom=0.155,
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
        ("The 1730 peak is ENGLISH.\n"
         "Same window, same method:  HadCET %+.2f K\n"
         "N. Hemisphere %+.3f K  (5-95%%: %+.3f to %+.3f)\n"
         "The English excursion is the %.0fth percentile of all %d\n"
         "55-year windows in the CET record. The hemisphere barely moves."
         % (exc_cet, exc_nh, d["exc_ci"]["nh"]["p5"], d["exc_ci"]["nh"]["p95"],
            d["cet_pct"], d["cet_windows_n"])),
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
    cy = d["compare_year"]
    names = ["CO$_2$\n(%d-%d)" % (FORCING_REF_YEAR, cy),
             "Solar\n(%d-%d)" % (FORCING_REF_YEAR, cy),
             "Solar, most\nfavourable case\n(%d-%d)" % (d["trough_year"], cy)]
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
    # HadCET is a single series and has no ensemble, so it carries no interval.
    # The reconstructions do, and hiding it would make a resolved and an
    # unresolved result look equally certain.
    lo = [0.0] + [d["exc"][k] - d["exc_ci"][k]["p5"] for k in ("nh", "sh")]
    hi = [0.0] + [d["exc_ci"][k]["p95"] - d["exc"][k] for k in ("nh", "sh")]
    bars = a3.bar([EXC_LABEL[k] for k in keys], vs,
                  color=[C_CET, C_REC, C_REC], width=0.6,
                  yerr=[lo, hi], capsize=5,
                  error_kw=dict(ecolor="#333333", lw=1.3))
    a3.axhline(0, color="#999999", lw=0.8, ls=":", zorder=0)
    for b_, k in zip(bars, keys):
        v = d["exc"][k]
        top = v if k == "cet" else d["exc_ci"][k]["p95"]
        a3.text(b_.get_x() + b_.get_width() / 2, top + 0.045, "%+.3f K" % v,
                ha="center", fontsize=9.0, fontweight="bold")
    a3.set_ylabel("Warming, 1680-1700 to 1725-1735 (K)", fontsize=9)
    a3.set_ylim(min(-0.16, d["exc_ci"]["sh"]["p5"] * 1.35), max(vs) * 1.45)
    a3.tick_params(labelsize=8)
    a3.set_title("The same event, measured the same way\n"
                 "England %+.2f K, hemispheres near zero"
                 % d["exc"]["cet"], fontsize=8.8, pad=6)
    a3.text(0.02, 0.975,
            "bars: median of %d members\nwhiskers: 5-95%% of members"
            % d["exc_ci"]["nh"]["n"],
            transform=a3.transAxes, ha="left", va="top", fontsize=6.9,
            color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", fc="#F2F6FA", ec="#9DB6CC"))

    for a_ in (ax, a0, a1, a2, a3):
        a_.grid(alpha=0.16, lw=0.6)
        a_.set_axisbelow(True)

    fig.text(0.055, 0.972, "Was it the Sun?", fontsize=21, fontweight="bold", va="top")
    fig.text(0.055, 0.934,
             "The same 1600-2025 chart as the version circulating on X, rebuilt with every "
             "series loaded from its published archive and every axis in comparable units.",
             fontsize=10.5, va="top", color="#333333")
    fig.text(0.055, 0.905,
             "Corrections: (1) raw TSI in W/m2 is not a radiative forcing - the "
             "conversion is dTSI/4x(1-albedo) = %.3f, so plotting raw dTSI "
             "overstates the solar term %.1fx on geometry alone.\n"
             "(2) Central England is not the globe.   "
             "(3) the Spoerer Minimum lies outside a 1600-start axis."
             % (TSI_TO_FORCING, 1.0 / TSI_TO_FORCING),
             fontsize=9, va="top", color="#555555")
    fig.text(0.945, 0.972,
             "Sources: HadCET · HadCRUT5 · UAH v6.1 · Law Dome (MacFarling Meure 2006) · "
             "NOAA Mauna Loa\nNRLTSI2 (Coddington 2016) · Quelccaya (Thompson 1986) · "
             "Neukom 2018 / PAGES 2k v2\nFull citations, URLs and SHA-256 of every file: "
             "MANIFEST in climate_figure.py",
             fontsize=7.8, va="top", ha="right", color="#444444", linespacing=1.6)

    fig.text(0.055, 0.100,
             "Hemispheric excursions are ensemble medians of %d members; whiskers "
             "and the ranges quoted above are the 5-95%% spread. %.0f%% of NH members "
             "are positive, but only %.0f%% of SH members are, so the SIGN of the "
             "Southern Hemisphere excursion is not established and is not claimed "
             "here.\nAcross %d alternative window definitions (start 1670-1690, baseline 15/20/25 yr, peak +40-50/+45-55/+50-60) the "
             "hemispheric excursion never exceeds %+.2f K (median %+.3f K), while the English one runs %+.2f to %+.2f K (median %+.2f K):\n"
             "the published window is at the favourable end for England; the hemispheric result is insensitive to the choice.\n"
             "The ratio of the two excursions is deliberately not quoted as "
             "the headline: its denominator is a hemispheric mean close to zero, "
             "which makes the quotient unstable. The two absolute numbers, and where "
             "the English one sits among all CET windows, are the stable statements."
             % (d["exc_ci"]["nh"]["n"], d["exc_ci"]["nh"]["frac_pos"] * 100,
                d["exc_ci"]["sh"]["frac_pos"] * 100,
                d["alt_n"], d["alt_nh"][2], d["alt_nh"][0],
                d["alt_cet"][1], d["alt_cet"][2], d["alt_cet"][0]),
             fontsize=7.8, va="top", color="#444444", linespacing=1.5)
    fig.text(0.055, 0.034,
             "Quelccaya ice cap, Peru is the source figure's own reference (3) for "
             "\"seen in many locations globally\". Over this window its d18O changes by "
             "%+.2f permil. Tropical d18O is precipitation-influenced,\nso it is not read "
             "here as a temperature in either direction: the record simply provides no "
             "clear evidence of a parallel global warming."
             % d["exc_quel"],
             fontsize=7.8, va="top", color="#444444", linespacing=1.5)

    fig.savefig(out, dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
