"""Download every raw file the figure needs into data/.

Nothing here transforms data. It fetches the exact files listed in the MANIFEST
of climate_figure.py, reports the sha256 of what arrived, and compares it with
the hash the pipeline expects.

A hash mismatch is NOT an error here -- providers do reissue datasets. It is a
signal to go and find out what changed before touching EXPECT. This script will
never edit EXPECT for you.

    python fetch_data.py            # fetch anything missing
    python fetch_data.py --force    # re-fetch everything
"""

from __future__ import annotations

import argparse
import sys
import urllib.error
import urllib.request
from pathlib import Path

from climate_figure import CFG, DATA, EXPECT, sha256_of

# Exact file URLs, matching the MANIFEST entries one for one.
URLS = {
    "hadcet": "https://www.metoffice.gov.uk/hadobs/hadcet/data/"
              "meantemp_monthly_totals.txt",
    "hadcrut5": "https://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.0.2.0/"
                "analysis/diagnostics/"
                "HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv",
    "uah_lt": "https://www.nsstc.uah.edu/data/msu/v6.1/tlt/uahncdc_lt_6.1.txt",
    "co2_ice": "https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/law/"
               "law2006.txt",
    "co2_mlo": "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_mlo.txt",
    "tsi": "https://lasp.colorado.edu/lisird/latis/dap/nrl2_tsi_P1Y.csv",
    "quelccaya": "https://www.ncei.noaa.gov/pub/data/paleo/icecore/trop/quelccaya/"
                 "q83summ-noaa.txt",
    "recon_nh": "https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/neukom2018/"
                "Real_proxy_recons/NH.txt",
    "recon_sh": "https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/neukom2018/"
                "Real_proxy_recons/SH.txt",
    "sunspots": "https://www.sidc.be/SILSO/DATA/SN_y_tot_V2.0.txt",
}

# GloSAT is deliberately absent: it needs a free CEDA account and the figure uses
# HadCRUT5 instead, labelled as HadCRUT5. Do not substitute one for the other.


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="re-download files that are already present")
    ap.add_argument("--timeout", type=int, default=120)
    args = ap.parse_args()

    missing = sorted(set(CFG) - set(URLS))
    if missing:
        print(f"No URL for: {', '.join(missing)}. Add it here and to the MANIFEST.",
              file=sys.stderr)
        return 1

    DATA.mkdir(exist_ok=True)
    ok = changed = failed = 0

    for key in CFG:
        dest = DATA / CFG[key]["file"]
        if dest.exists() and not args.force:
            print(f"  have    {key:<10s} {dest.name}")
            ok += 1
            continue
        try:
            print(f"  fetch   {key:<10s} {URLS[key]}")
            with urllib.request.urlopen(URLS[key], timeout=args.timeout) as r:
                blob = r.read()
            dest.write_bytes(blob)
        except (urllib.error.URLError, OSError) as exc:
            print(f"  FAILED  {key:<10s} {exc}", file=sys.stderr)
            failed += 1
            continue

        got = sha256_of(dest)
        want = EXPECT[key]["sha256"]
        if got == want:
            print(f"          sha256 matches EXPECT ({got[:12]})")
            ok += 1
        else:
            print(f"          sha256 CHANGED\n"
                  f"            expected {want}\n"
                  f"            got      {got}\n"
                  f"          The provider reissued this dataset. Find out what "
                  f"changed before updating EXPECT.")
            changed += 1

    print(f"\n{ok} ok, {changed} changed, {failed} failed.")
    if failed:
        print("A missing file is a stop, not a licence to substitute another one.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
