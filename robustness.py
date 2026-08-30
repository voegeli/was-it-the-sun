"""Sensitivity of the CO2:solar forcing ratio to every parameter an opponent
can plausibly contest. Run alongside climate_figure.py; prints a table.

The point is not that our numbers are exact. It is that the conclusion does not
depend on the choices we made -- it survives every reasonable alternative.
"""

import numpy as np

import climate_figure as cf

series, _, _ = cf.load_all()
d = cf.derive(series)

tsi = d["tsi"]
co2 = d["co2"]
# Both terms are read at the same year -- the last one where the CO2 record and
# the centred 11-year solar mean are both defined. Reading CO2 at 2025 against a
# solar mean that ends 2018 would inflate every ratio below.
CY = d["compare_year"]
c_now = float(co2.loc[CY])
t_now = float(tsi.loc[CY])
print(f"\nAll comparisons evaluated at the common end year {CY} "
      f"(CO2 data to {int(co2.index.max())}, centred 11-yr solar mean to "
      f"{int(tsi.index.max())}).")

print("\n" + "=" * 78)
print("ROBUSTNESS OF THE CO2 : SOLAR FORCING RATIO")
print("=" * 78)

print("\n1. Albedo (we use 0.30). Factor = (1-alpha)/4")
for alb in (0.00, 0.15, 0.30, 0.35):
    f = (1 - alb) / 4
    fs = (t_now - float(tsi.loc[1750])) * f
    fc = cf.CO2_FORCING_COEFF * np.log(c_now / float(series["co2_ice"].loc[1750]))
    print(f"   alpha={alb:.2f}  factor={f:.4f}  solar={fs:+.4f}  ratio={fc/fs:6.1f} : 1")

print("\n2. Reference year for BOTH terms (we use 1750, per IPCC AR6)")
for ref in (1620, 1650, 1700, 1750, 1800, 1850):  # 11-yr mean starts 1615
    fs = (t_now - float(tsi.loc[ref])) * cf.TSI_TO_FORCING
    fc = cf.CO2_FORCING_COEFF * np.log(c_now / float(series["co2_ice"].loc[ref]))
    print(f"   ref={ref}  solar={fs:+.4f}  CO2={fc:+.4f}  ratio={fc/fs:6.1f} : 1")

print("\n3. Solar measured from its single lowest point (maximally favourable)")
ty = int(tsi.loc[1600:1720].idxmin())
fs = (t_now - float(tsi.loc[ty])) * cf.TSI_TO_FORCING
fc = cf.CO2_FORCING_COEFF * np.log(c_now / float(series["co2_ice"].loc[1750]))
print(f"   solar from {ty} vs CO2 from 1750, both to {CY}: ratio={fc/fs:6.1f} : 1")

print("\n4. Raw TSI difference, NO conversion at all (the source figure's method)")
raw = t_now - float(tsi.loc[1750])
print(f"   raw dTSI={raw:+.4f} W/m2 vs CO2 forcing {fc:+.4f} W/m2  ratio={fc/raw:6.1f} : 1")
print("   -- even skipping the geometry entirely, CO2 still wins.")

print("\n5. Cross-check of our CO2 forcing against IPCC AR6 WG1 Ch.7")
fc19 = cf.CO2_FORCING_COEFF * np.log(
    float(co2.loc[2019]) / float(series["co2_ice"].loc[1750]))
print(f"   ours, 1750-2019: {fc19:+.3f} W/m2   AR6 published: +2.16 W/m2")

print("\n6. The sign test -- needs no amplitude argument at all")
print(f"   solar forcing trend since 1980: {d['tr_fsol'][0]:+.4f} W/m2/decade "
      f"({d['tr_fsol'][1]}-{d['tr_fsol'][2]})")
print(f"   HadCRUT5 trend since 1980:      {d['tr_had'][0]:+.4f} K/decade")
print(f"   UAH trend since 1980:           {d['tr_uah'][0]:+.4f} K/decade")
print("   Opposite signs. No rescaling of either axis can repair a sign.")

print("\n7. Amplification needed to close the gap")
need = d["f_co2_now"] / d["f_sol_now"]
print(f"   a solar amplifier would have to multiply the solar term by {need:.0f}x")
print("   AND reverse the sign of its post-1980 trend to explain the warming.")
print("=" * 78 + "\n")

print("\n8. THE RECONSTRUCTION-CHOICE ATTACK: what if NRLTSI2 is too flat?")
# The predictable objection is that this figure picked a low-variability TSI
# reconstruction. Rather than argue about which reconstruction is right, compute
# what that objection would have to be worth to change the conclusion.
dtsi_1750 = t_now - float(tsi.loc[1750])
need_factor = fc / (dtsi_1750 * cf.TSI_TO_FORCING)
need_dtsi = dtsi_1750 * need_factor
print(f"   NRLTSI2 dTSI 1750-{CY}: {dtsi_1750:+.4f} W/m2")
print("   For the solar forcing to EQUAL the CO2 forcing, dTSI would have to be")
print(f"   {need_factor:.0f}x larger = {need_dtsi:.1f} W/m2, i.e. "
      f"{need_dtsi / t_now * 100:.2f}% of the solar constant.")
print("   A critic who rejects NRLTSI2 has to name a reconstruction that reaches")
print("   that, not merely one that is somewhat higher.")
print("\n   The sign test is LESS reconstruction-dependent than the amplitude")
print("   comparison, but it is not independent. NRLTSI2 is calibrated against")
print("   direct space-based observations over the satellite era and reproduces")
print("   their variability, so its modern trend is far better constrained than")
print("   its Maunder-era amplitude -- but it is still a proxy model built on")
print("   sunspot and facular indices, and no alternative TSI composite is loaded")
print("   here. The claim this pipeline can support is narrow: within NRLTSI2 the")
print("   trend since 1980 is negative.")
print("\n" + "=" * 78 + "\n")
