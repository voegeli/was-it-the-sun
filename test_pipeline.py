"""Targeted tests for the guarantees this pipeline claims to make.

These are not unit tests of plotting. They test the four things a reader of the
figure has to be able to trust: that the two forcings are compared at the same
instant, that the reconstruction uncertainty is real and carried through, that
the MANIFEST is provenance rather than prose, and that no series is loaded and
then quietly ignored.

Runs under pytest, or standalone:  python test_pipeline.py
"""

from __future__ import annotations

import copy
import datetime

import numpy as np

import climate_figure as cf

_SERIES = None
_DERIVED = None


def loaded():
    global _SERIES, _DERIVED
    if _SERIES is None:
        _SERIES, _, _ = cf.load_all()
        _DERIVED = cf.derive(_SERIES)
    return _SERIES, _DERIVED


# --------------------------------------------------------------------------
# 1. common end year for the forcing comparison
# --------------------------------------------------------------------------


def test_compare_year_is_defined_in_both_series():
    _, d = loaded()
    cy = d["compare_year"]
    assert cy in d["f_co2"].index, "compare_year missing from the CO2 forcing series"
    assert cy in d["f_sol"].index, "compare_year missing from the solar forcing series"


def test_compare_year_is_the_last_year_both_are_defined():
    _, d = loaded()
    expected = int(min(d["f_co2"].index.max(), d["f_sol"].index.max()))
    assert d["compare_year"] == expected


def test_ratio_uses_the_same_year_for_both_terms():
    _, d = loaded()
    cy = d["compare_year"]
    assert d["f_co2_now"] == float(d["f_co2"].loc[cy])
    assert d["f_sol_now"] == float(d["f_sol"].loc[cy])
    # the bug this guards: CO2 read at the end of its record, solar at the end of
    # the shorter centred mean, both labelled "today"
    assert d["f_co2_now"] != float(d["f_co2"].iloc[-1]) or cy == int(
        d["f_co2"].index.max()
    )


def test_solar_mean_is_centred_not_trailing():
    """A trailing mean would silently extend the solar series to the CO2 end year."""
    _, d = loaded()
    raw_end = int(d["tsi_raw"].index.max())
    smooth_end = int(d["f_sol"].index.max())
    assert smooth_end < raw_end, "centred 11-yr mean must stop short of the raw record"
    assert raw_end - smooth_end == 5, "an 11-year centred window loses 5 years at each end"


# --------------------------------------------------------------------------
# 2. ensemble uncertainty
# --------------------------------------------------------------------------


def test_ensemble_intervals_exist_and_are_ordered():
    _, d = loaded()
    for key in ("nh", "sh"):
        ci = d["exc_ci"][key]
        assert ci["n"] == 100, f"{key}: expected 100 ensemble members, got {ci['n']}"
        assert ci["p5"] <= ci["median"] <= ci["p95"], f"{key}: percentiles out of order"
        assert 0.0 <= ci["frac_pos"] <= 1.0


def test_plotted_bar_and_interval_are_the_same_statistic():
    """The bar must be the median OF the excursions, not the excursion OF the median."""
    _, d = loaded()
    for key in ("nh", "sh"):
        assert d["exc"][key] == d["exc_ci"][key]["median"]


def test_intervals_are_recomputed_from_the_ensemble_files():
    _, d = loaded()
    for key in ("nh", "sh"):
        path, _ = cf.resolve("recon_" + key)
        ens = cf._load_recon_ensemble(path)
        diff = (ens.loc[1725:1735].mean(axis=0) - ens.loc[1680:1700].mean(axis=0))
        diff = diff.to_numpy(dtype=float)
        ci = d["exc_ci"][key]
        assert np.isclose(ci["median"], float(np.median(diff)))
        assert np.isclose(ci["p5"], float(np.percentile(diff, 5)))
        assert np.isclose(ci["p95"], float(np.percentile(diff, 95)))


def test_southern_hemisphere_sign_is_not_claimed():
    """SH interval spans zero, so nothing downstream may assert its sign."""
    _, d = loaded()
    ci = d["exc_ci"]["sh"]
    if ci["p5"] < 0 < ci["p95"]:
        assert ci["frac_pos"] < 0.95, (
            "SH interval spans zero but frac_pos suggests otherwise; check the maths"
        )


def test_ratio_is_computed_but_not_the_headline():
    """The quotient is unstable -- its denominator is near zero in many windows.

    It stays in the QC output for anyone who wants it, but must not appear in a
    panel title or in the main callout, where it would read as the finding.
    """
    _, d = loaded()
    assert d["exc_ratio_min"] == d["exc"]["cet"] / d["exc_ci"]["nh"]["p95"]
    src = cf.Path(cf.__file__).read_text(encoding="utf-8")
    fig = src[src.index("def make_figure"):]
    for banned in ['exc_ratio_min"]', 'exc_ratio"]']:
        assert banned not in fig, (
            "the excursion ratio leaked back into the figure; lead with the two "
            "absolute numbers instead"
        )


def test_cet_percentile_is_computed_from_all_windows():
    _, d = loaded()
    cet = d["cet"]
    win = []
    for y0 in range(int(cet.index.min()), int(cet.index.max()) - 55):
        a = cet.loc[y0:y0 + 20].mean()
        b = cet.loc[y0 + 45:y0 + 55].mean()
        if np.isfinite(a) and np.isfinite(b):
            win.append(b - a)
    win = np.asarray(win, dtype=float)
    assert d["cet_windows_n"] == win.size
    assert np.isclose(d["cet_pct"], (win < d["exc"]["cet"]).mean() * 100.0)
    assert 0.0 <= d["cet_pct"] <= 100.0


def test_audit_and_pipeline_share_one_window_population():
    """audit.py must measure the same windows the published figure quotes."""
    _, d = loaded()
    win = cf.excursion_windows(d["cet"])
    assert len(win) == d["cet_windows_n"]
    assert int(win.index.min()) == d["cet_windows_first"]
    assert int(win.index.max()) == d["cet_windows_last"]
    pct = float((win.to_numpy() < d["exc"]["cet"]).mean() * 100.0)
    assert np.isclose(pct, d["cet_pct"], atol=1e-12)


def test_window_population_starts_at_the_first_year_of_the_record():
    """No hand-picked lower bound: that is how the two counts diverged before."""
    _, d = loaded()
    assert d["cet_windows_first"] == int(d["cet"].index.min())
    assert d["cet_windows_last"] == int(d["cet"].index.max()) - cf.EXC_LEN - 1


def test_published_window_is_a_member_of_the_family():
    """1680-1700 vs 1725-1735 must be exactly the y0 = EXC_START window."""
    _, d = loaded()
    assert cf.EXC_BASELINE == (0, 20) and cf.EXC_PEAK == (45, 55)
    assert cf.EXC_START + cf.EXC_BASELINE[1] == 1700
    assert cf.EXC_START + cf.EXC_PEAK[0] == 1725
    assert cf.EXC_START + cf.EXC_PEAK[1] == 1735
    assert np.isclose(cf.excursion(d["cet"], cf.EXC_START), d["exc"]["cet"])


def test_audit_has_no_hardcoded_window_start():
    src = cf.Path("audit.py").read_text(encoding="utf-8")
    block = src[src.index("THE ABSOLUTE CLAIM"):src.index("CO2 REFERENCE VALUE")]
    assert "1660" not in block, "audit still hardcodes a window start year"
    assert "excursion_windows" in block


def test_co2_reference_is_loaded_not_a_literal():
    """C0 must come from the ice core, not from a typed-in pre-industrial value."""
    series, d = loaded()
    c0 = float(series["co2_ice"].loc[cf.FORCING_REF_YEAR])
    expected = cf.CO2_FORCING_COEFF * np.log(
        float(d["co2"].loc[d["compare_year"]]) / c0)
    assert np.isclose(d["f_co2_now"], expected)
    src = cf.Path(cf.__file__).read_text(encoding="utf-8")
    assert "278.3" not in src, "a published pre-industrial value was hardcoded"


# --------------------------------------------------------------------------
# README consistency
# --------------------------------------------------------------------------


def test_readme_headline_numbers_match_the_run():
    """The README transcribes rounded values by hand. Catch it when they drift.

    Deliberately not a README-rewriting system: this only asserts that what is
    written there still matches what derive() produces, at the precision the
    README quotes.
    """
    _, d = loaded()
    readme = cf.Path("README.md").read_text(encoding="utf-8")
    expected = [
        (f"{cf.TSI_TO_FORCING:.3f}", "TSI-to-forcing factor"),
        (f"{1.0 / cf.TSI_TO_FORCING:.1f}", "geometric overstatement factor"),
        (f"{d['compare_year']}", "common end year"),
        (f"{d['f_co2_now']:+.3f}".replace("+", ""), "CO2 forcing"),
        (f"{d['f_sol_now']:+.3f}".replace("+", ""), "solar forcing"),
        (f"{d['ratio']:.0f} : 1", "forcing ratio"),
        (f"{d['ratio_best']:.0f} : 1", "best-case ratio"),
        (f"{d['exc']['cet']:+.3f}".replace("+", ""), "HadCET excursion"),
        (f"{d['exc']['nh']:+.3f}".replace("+", ""), "NH excursion"),
        (f"{d['exc']['sh']:+.3f}".replace("+", ""), "SH excursion"),
        (f"{d['cet_windows_n']}", "CET window count"),
        (f"{d['cet_pct']:.0f}th percentile", "CET percentile"),
    ]
    missing = [(v, what) for v, what in expected if v not in readme]
    assert not missing, (
        "README no longer matches the run: " +
        ", ".join(f"{what}={v!r}" for v, what in missing)
    )


def test_readme_does_not_claim_its_own_numbers_are_computed():
    readme = cf.Path("README.md").read_text(encoding="utf-8")
    assert "none are typed in" not in readme, (
        "the README quotes hand-transcribed values; it must not claim otherwise"
    )


# --------------------------------------------------------------------------
# 3. MANIFEST is machine-checkable provenance
# --------------------------------------------------------------------------


def test_manifest_cfg_expect_cover_the_same_series():
    assert set(cf.MANIFEST) == set(cf.CFG) == set(cf.EXPECT)


def test_every_manifest_field_present_and_non_empty():
    for key, entry in cf.MANIFEST.items():
        for field in cf.REQUIRED_MANIFEST_FIELDS:
            assert field in entry, f"{key}: missing {field}"
            assert isinstance(entry[field], str) and entry[field].strip(), (
                f"{key}: {field} is empty"
            )


def test_manifest_urls_are_urls_and_file_url_is_not_the_landing_page():
    for key, entry in cf.MANIFEST.items():
        for field in ("landing", "file_url"):
            assert entry[field].startswith(("http://", "https://")), f"{key}: {field}"
        assert entry["file_url"].rstrip("/") != entry["landing"].rstrip("/"), key


def test_access_dates_parse():
    for key, entry in cf.MANIFEST.items():
        datetime.date.fromisoformat(entry["accessed"])


def test_sha256_lives_in_the_manifest_and_only_there():
    """One hash, one home. Two hand-maintained copies drift."""
    for key, entry in cf.MANIFEST.items():
        assert len(entry["sha256"]) == 64, f"{key}: sha256 is not 64 hex characters"
        int(entry["sha256"], 16)
        assert entry["sha256"] == entry["sha256"].lower()
    for key, exp in cf.EXPECT.items():
        assert "sha256" not in exp, (
            f"EXPECT[{key!r}] still carries a duplicate sha256"
        )


def _manifest_mutation_raises(key, field, value, delete=False):
    original = copy.deepcopy(cf.MANIFEST)
    try:
        if delete:
            del cf.MANIFEST[key][field]
        else:
            cf.MANIFEST[key][field] = value
        try:
            cf.check_manifest_completeness()
        except cf.DataError:
            return True
        return False
    finally:
        cf.MANIFEST.clear()
        cf.MANIFEST.update(original)


def test_missing_sha256_is_rejected():
    assert _manifest_mutation_raises("tsi", "sha256", None, delete=True)


def test_empty_sha256_is_rejected():
    assert _manifest_mutation_raises("tsi", "sha256", "   ")


def test_short_sha256_is_rejected():
    assert _manifest_mutation_raises("tsi", "sha256", "abc123")


def test_long_sha256_is_rejected():
    assert _manifest_mutation_raises("tsi", "sha256", "a" * 65)


def test_non_hex_sha256_is_rejected():
    assert _manifest_mutation_raises("tsi", "sha256", "z" * 64)


def test_uppercase_sha256_is_rejected():
    assert _manifest_mutation_raises(
        "tsi", "sha256", cf.MANIFEST["tsi"]["sha256"].upper())


def test_validation_compares_against_the_manifest_hash():
    """Flipping the MANIFEST hash must make loading fail, proving it is the one used."""
    original = copy.deepcopy(cf.MANIFEST)
    try:
        cf.MANIFEST["tsi"]["sha256"] = "0" * 64
        raised = False
        try:
            cf.load_all()
        except cf.DataError as exc:
            raised = "sha256 mismatch" in str(exc)
        assert raised, "validate() is not checking against MANIFEST[key]['sha256']"
    finally:
        cf.MANIFEST.clear()
        cf.MANIFEST.update(original)


def test_downloader_checks_the_manifest_hash():
    import fetch_data
    src = cf.Path(fetch_data.__file__).read_text(encoding="utf-8")
    assert 'MANIFEST[key]["sha256"]' in src
    assert "EXPECT" not in src, "fetch_data still reads hashes from EXPECT"


def test_empty_citation_is_rejected():
    original = copy.deepcopy(cf.MANIFEST)
    try:
        cf.MANIFEST["tsi"]["citation"] = "   "
        raised = False
        try:
            cf.check_manifest_completeness()
        except cf.DataError:
            raised = True
        assert raised, "a blank citation must raise DataError"
    finally:
        cf.MANIFEST.clear()
        cf.MANIFEST.update(original)


def test_missing_url_is_rejected():
    original = copy.deepcopy(cf.MANIFEST)
    try:
        del cf.MANIFEST["tsi"]["file_url"]
        raised = False
        try:
            cf.check_manifest_completeness()
        except cf.DataError:
            raised = True
        assert raised, "a missing file_url must raise DataError"
    finally:
        cf.MANIFEST.clear()
        cf.MANIFEST.update(original)


def test_orphan_manifest_entry_is_rejected():
    original = copy.deepcopy(cf.MANIFEST)
    try:
        cf.MANIFEST["not_a_series"] = dict(cf.MANIFEST["tsi"])
        raised = False
        try:
            cf.check_manifest_completeness()
        except cf.DataError:
            raised = True
        assert raised, "a MANIFEST entry with no CFG series must raise DataError"
    finally:
        cf.MANIFEST.clear()
        cf.MANIFEST.update(original)


def test_downloader_uses_the_manifest_urls():
    """The onaa.gov/noaa.gov class of bug: two copies of one URL, one of them wrong."""
    import fetch_data

    assert set(fetch_data.URLS) == set(cf.CFG)
    for key in cf.CFG:
        assert fetch_data.URLS[key] == cf.MANIFEST[key]["file_url"]


def test_declared_columns_match_the_real_headers():
    for key in cf.CFG:
        path, _ = cf.resolve(key)
        cf.check_column(key, path)  # raises DataError on drift


# --------------------------------------------------------------------------
# 4./5. no half-registered series
# --------------------------------------------------------------------------


def test_sunspots_fully_removed():
    for container in (cf.CFG, cf.EXPECT, cf.MANIFEST, cf.LOADERS):
        assert "sunspots" not in container
    assert not hasattr(cf, "load_sunspots")


def test_every_loaded_series_is_used_downstream():
    series, d = loaded()
    used = {
        "hadcet": "cet", "hadcrut5": "hadcrut5", "uah_lt": "uah",
        "co2_ice": "f_co2_ice", "co2_mlo": "f_co2", "tsi": "f_sol",
        "quelccaya": "exc_quel", "recon_nh": "nh", "recon_sh": "sh",
    }
    assert set(used) == set(series), (
        "a series is loaded and validated but never used, or vice versa"
    )
    for key in used.values():
        assert key in d, f"derived value {key} missing"


def test_quelccaya_is_never_converted_to_a_temperature():
    """Its sign must not be reported as warming or cooling anywhere in the source."""
    src = (cf.Path(cf.__file__)).read_text(encoding="utf-8")
    banned = ["opposite direction", "Quelccaya cooling", "cooling in Peru"]
    for phrase in banned:
        assert phrase not in src, f"temperature reading of d18O leaked back in: {phrase}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except AssertionError as exc:
            failed += 1
            print(f"  FAIL  {fn.__name__}: {exc}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"  ERROR {fn.__name__}: {type(exc).__name__}: {exc}")
    print(f"\n{len(fns) - failed} passed, {failed} failed")
    raise SystemExit(1 if failed else 0)
