"""Download every raw file the figure needs into data/.

Nothing here transforms data. It fetches the exact files listed in the MANIFEST
of climate_figure.py, reports the sha256 of what arrived, and compares it with
the hash the pipeline expects.

A hash mismatch is NOT an error here -- providers do reissue datasets. It is a
signal to go and find out what changed before touching MANIFEST. This script
will never edit MANIFEST for you.

    python fetch_data.py            # fetch anything missing
    python fetch_data.py --force    # re-fetch everything
"""

from __future__ import annotations

import argparse
import sys
import urllib.error
import urllib.request

from climate_figure import CFG, DATA, MANIFEST, sha256_of

# URLs are NOT duplicated here. They come from MANIFEST in climate_figure.py, so
# the documented source and the fetched source cannot drift apart -- an earlier
# revision kept a private copy of this list and one entry silently disagreed with
# the MANIFEST by a single character in the hostname.
URLS = {key: MANIFEST[key]["file_url"] for key in CFG}

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
        want = MANIFEST[key]["sha256"]
        if got == want:
            print(f"          sha256 matches MANIFEST ({got[:12]})")
            ok += 1
        else:
            print(f"          sha256 CHANGED\n"
                  f"            expected {want}\n"
                  f"            got      {got}\n"
                  f"          The provider reissued this dataset. Find out what "
                  f"changed before updating MANIFEST.")
            changed += 1

    print(f"\n{ok} ok, {changed} changed, {failed} failed.")
    if failed:
        print("A missing file is a stop, not a licence to substitute another one.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
