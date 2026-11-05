
PLL Loop-filter design Script - 2.4 GHz Bluetooth LO synthesizer

Specs in = loop-filter Component values out, plus a check of the
achieved Bandwidth and phase margin using the exact filter impedance.

      (3rd-order passive filter, 4th-order type-II loop):

   charge pump ──┬────────┬──── R3 ────┬──── VCO tune (Vc)
                 C1       R2           C3
                 │        C2           │
                GND      GND          GND

Run Command:  python3 pll_design.py
"""
import math, cmath

#SPECS
FREF  = 1e6          # channel spacing (Hz) 
N_MIN, N_MAX = 2402, 2480   # BLE/BT channels 2402-2480 (MHz)
F0    = 2.1e9        # VCO free-running freq at Vc = 0 (Hz)
KV    = 150e6        # VCO gain (Hz/V)
ICP   = 1e-3         # charge-pump current (A)
FC    = 100e3        # loop bandwidth (Hz)  
PM    = 55.0         # phase margin (deg)    
T31   = 0.25         # 3rd pole time constant/fraction of T1 
C3    = 10e-12       # chosen; R3 = T3/C3 


N   = round((N_MIN + N_MAX) / 2)    
wc  = 2 * math.pi * FC
Kphi = ICP / (2 * math.pi)        
Kvco = 2 * math.pi * KV              
pm   = math.radians(PM)

def phase_margin(T1):
    T3 = T31 * T1
    T2 = 1 / (wc**2 * (T1 + T3))    
    return math.atan(wc*T2) - math.atan(wc*T1) - math.atan(wc*T3), T2, T3

# bisection for T1 that gives the requested phase margin
lo, hi = 1e-12, 1 / wc
for _ in range(200):
    mid = (lo + hi) / 2
    if phase_margin(mid)[0] > pm: lo = mid
    else: hi = mid
T1 = mid
_, T2, T3 = phase_margin(T1)

Ctot = (Kphi * Kvco / (N * wc**2)) * math.sqrt(1 + (wc*T2)**2) / (
        math.sqrt(1 + (wc*T1)**2) * math.sqrt(1 + (wc*T3)**2))
C1 = Ctot * T1 / T2
C2 = Ctot - C1
R2 = T2 / C2
R3 = T3 / C3

#  exact check with the real impedance including C3/R3 loading
def Zfilt(s):
    zc1 = 1/(s*C1); zb = R2 + 1/(s*C2); zc3 = 1/(s*C3)
    z_in = 1/(1/zc1 + 1/zb + 1/(R3 + zc3))
    return z_in * zc3 / (R3 + zc3)           

def G(f, n=N):
    s = 2j*math.pi*f
    return Kphi * Zfilt(s) * Kvco / (s * n)

f = 1e3; prev = abs(G(f))
while abs(G(f)) > 1: f *= 1.0005
fc_real = f
pm_real = 180 + math.degrees(cmath.phase(G(fc_real)))
# reference-spur suppression by the extra R3C3 pole at FREF
pole3 = 1/(2*math.pi*R3*C3)
extra = 20*math.log10(abs(1/(1+1j*FREF/pole3)))

print("=== PLL loop-filter design (2.4 GHz BLE synthesizer) ===")
print(f"N (design) = {N},  Kphi = {Kphi*1e6:.1f} uA/rad,  Kvco = 2pi x {KV/1e6:.0f} MHz/V")
print(f"T1 = {T1*1e9:.1f} ns  T2 = {T2*1e6:.3f} us  T3 = {T3*1e9:.1f} ns")
print(f"C1 = {C1*1e12:.1f} pF")
print(f"C2 = {C2*1e12:.1f} pF")
print(f"R2 = {R2/1e3:.2f} kOhm")
print(f"R3 = {R3/1e3:.2f} kOhm   C3 = {C3*1e12:.1f} pF  (extra pole at {pole3/1e6:.2f} MHz, {extra:.1f} dB more spur rejection at FREF)")
print(f"Exact check: unity-gain BW = {fc_real/1e3:.1f} kHz, phase margin = {pm_real:.1f} deg")
for n in (N_MIN, N_MAX):
    print(f"   at N={n}: Vc = {(n*FREF-F0)/KV:.3f} V")
# rounded to standard E-series values (used in the SPICE netlists)
C1, C2, R2, R3, C3 = 39e-12, 470e-12, 11e3, 10e3, 10e-12
f = 1e3
while abs(G(f)) > 1: f *= 1.0005
print(f"With E-series values C1=39p C2=470p R2=11k R3=10k C3=10p: BW = {f/1e3:.1f} kHz, PM = {180+math.degrees(cmath.phase(G(f))):.1f} deg")
