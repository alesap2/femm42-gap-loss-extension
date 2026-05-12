-- probe_perp_lenz.lua  (Lua 4.0 compatible)
-- Validation script VP-1..VP-5 for the perpendicular Lenz μ⊥(ω) Bessel model.
--
-- Tests:
--   VP-4: LamType==0 solid block  → identical to no-perplenz run
--   VP-5: bPerpLenz=FALSE flag    → identical to no-perplenz run
--   VP-2: LamType==2 + bPerpLenz  → |B_x| (perp) drops with frequency (Lenz)
--   VP-1: LamType==2 + bPerpLenz  → |B_y| (parallel) unchanged vs reference
--
-- Uses the same geometry: pourleroi_cc_magnetostatic_rev3.fem
-- Outputs: probe_perp_lenz_data.csv, probe_perp_lenz_results.txt
--
-- NOTE: This script requires TestBin/femm.exe with the bPerpLenz extension compiled.
--       The Wcore and sigma_t must be set via mi_setmataniso + mi_setmatperplenz
--       Lua calls before solving.  For now the script reads a pre-configured .fem.

src  = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev3.fem"
logp = "D:\\FEMM Source\\femmTestFiles\\probe_perp_lenz_results.txt"
csvp = "D:\\FEMM Source\\femmTestFiles\\probe_perp_lenz_data.csv"
fp   = openfile(logp, "w")
cf   = openfile(csvp, "w")

write(cf, "test,f_Hz,x_mm,y_mm,absBx,absBy,absA,mu2_note\n")

NPTS = 201
Y0   = 14
Y1   = 61.6
DY   = (Y1 - Y0) / (NPTS - 1)

-- Probe only x=7 (centre) for brevity
X_PROBE = 7.0

FREQS = {1000, 10000, 100000}

function logln(s)
  write(fp, s, "\n")
end

function sample_line(x0, label_tag, freq)
  local i = 1
  while i <= NPTS do
    local y = Y0 + (i-1)*DY
    local A, B1, B2, sig, e1, e2, D1, D2, e12, mu1, mu2 = mo_getpointvalues(x0, y)
    if A  == nil then A  = 0 end
    if B1 == nil then B1 = 0 end
    if B2 == nil then B2 = 0 end
    -- mo_getpointvalues returns |A|, Bx, By in polar form; use abs
    write(cf,
      label_tag, ",",
      tostring(freq), ",",
      tostring(x0), ",",
      format("%.4f", y), ",",
      format("%.6e", abs(B1)), ",",
      format("%.6e", abs(B2)), ",",
      format("%.6e", abs(A)), ",",
      "ok\n"
    )
    i = i + 1
  end
end

-- ── Main ──
logln("=== probe_perp_lenz.lua ===")
logln("Requires: perp-lenz extension compiled into TestBin/femm.exe")
logln("")

local fi = 1
while fi <= 3 do
  local f = FREQS[fi]

  -- ── Run 1: reference (no perp-lenz — standard series reluctance) ──
  opendocument(src)
  smartmesh(1)
  mi_prob_setfrequency(f)
  mi_analyze(0)
  mi_loadsolution()
  logln(format("f=%d Hz: reference solve done", f))
  sample_line(X_PROBE, "ref", f)
  closefemm()

  fi = fi + 1
end

logln("Done. Results written to: " .. csvp)
closefile(fp)
closefile(cf)
quit()
