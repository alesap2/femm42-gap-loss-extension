-- probe_a_compare.lua  (Lua 4.0 compatible)
-- Tests user's physical critique: does A_z change when we kill conductivity?
-- If A_z is essentially identical between (A) full anisotropic sigma and
-- (B) sigma=0 everywhere, then the current model has NO Lenz back-reaction
-- for the perpendicular flux component.
-- Probe: vertical lines x={0.1,7,13.9}, y in [14, 61.6] mm in rev3.fem.

src  = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev3.fem"
logp = "D:\\FEMM Source\\femmTestFiles\\probe_a_compare_results.txt"
csvp = "D:\\FEMM Source\\femmTestFiles\\probe_a_compare_data.csv"
fp   = openfile(logp, "w")
cf   = openfile(csvp, "w")
write(cf, "f_Hz,x_mm,y_mm,absA_full,absA_sig0,absBx_full,absBx_sig0,absBy_full,absBy_sig0\n")

NPTS = 401
Y0   = 14
Y1   = 61.6
X_LIST = {0.1, 7, 13.9}
NX = 3
DY   = (Y1 - Y0) / (NPTS - 1)

FREQS = {1000, 10000, 100000}

function logln(s)
  write(fp, s, "\n")
end

function sample_line(x0)
  local arr = {}
  local i = 1
  while i <= NPTS do
    local y = Y0 + (i-1)*DY
    local A, B1, B2 = mo_getpointvalues(x0, y)
    if A  == nil then A  = 0 end
    if B1 == nil then B1 = 0 end
    if B2 == nil then B2 = 0 end
    arr[i] = {x=x0, y=y, A=A, B1=B1, B2=B2}
    i = i + 1
  end
  return arr
end

function compare(arrA, arrB, label)
  local maxdA   = 0
  local rmsdA   = 0
  local maxdB   = 0
  local rmsdB   = 0
  local maxAref = 0
  local maxBref = 0
  local i = 1
  while i <= NPTS do
    local dA   = abs(arrA[i].A - arrB[i].A)
    -- For complex phasors |B|^2 = |Bx|^2 + |By|^2  (not Bx^2 + By^2)
    local aBxR = abs(arrA[i].B1)
    local aByR = abs(arrA[i].B2)
    local aBxT = abs(arrB[i].B1)
    local aByT = abs(arrB[i].B2)
    local Bref = sqrt(aBxR*aBxR + aByR*aByR)
    local Btst = sqrt(aBxT*aBxT + aByT*aByT)
    local dB   = abs(Bref - Btst)
    if dA > maxdA then maxdA = dA end
    if dB > maxdB then maxdB = dB end
    rmsdA = rmsdA + dA*dA
    rmsdB = rmsdB + dB*dB
    local Aabs = abs(arrA[i].A)
    if Aabs > maxAref then maxAref = Aabs end
    if Bref > maxBref then maxBref = Bref end
    i = i + 1
  end
  rmsdA = sqrt(rmsdA / NPTS)
  rmsdB = sqrt(rmsdB / NPTS)
  local pctA = 0
  local pctB = 0
  if maxAref > 0 then pctA = 100*maxdA/maxAref end
  if maxBref > 0 then pctB = 100*maxdB/maxBref end
  logln(format("  %s :", label))
  logln(format("    max|A_full|                = %.6e Wb/m", maxAref))
  logln(format("    max|B_full|                = %.6e T",   maxBref))
  logln(format("    max|A_full - A_sigma=0|    = %.6e (%.3f %% of max|A|)", maxdA, pctA))
  logln(format("    rms|A_full - A_sigma=0|    = %.6e", rmsdA))
  logln(format("    max||B_full|-|B_sigma=0||  = %.6e (%.3f %% of max|B|)", maxdB, pctB))
  logln(format("    rms||B_full|-|B_sigma=0||  = %.6e", rmsdB))
end

logln("=== Probe A_z along x={0.1,7,13.9}, y in [14, 61.6] mm ===")
logln(format("Source: %s", src))
logln(format("Samples: %d points, dy = %.3f mm", NPTS, DY))
logln("Probe lines: x = 0.1, 7, 13.9 mm")
logln("Test: A_z(full anisotropic sigma) vs A_z(sigma=0 everywhere).")
logln("If results are nearly identical, the model has NO Lenz back-reaction:")
logln("losses are computed post-hoc from the macroscopic B that exists in a")
logln("pure magnetostatic-AC solve.")
logln("")

fi = 1
while fi <= 3 do
  f = FREQS[fi]
  logln(format("---- Frequency = %d Hz ----", f))

  -- Pass A: as-is (anisotropic sigma_t in gaps and core)
  open(src)
  mi_saveas("D:\\FEMM Source\\femmTestFiles\\probe_passA.fem")
  mi_probdef(f, "millimeters", "planar", 1e-008, 35, 30)
  mi_analyse(1)
  mi_loadsolution()
  arrA_all = {}
  xi = 1
  while xi <= NX do
    arrA_all[xi] = sample_line(X_LIST[xi])
    xi = xi + 1
  end
  mo_close()
  mi_close()

  -- Pass B: kill all conductivity in core and gap
  open(src)
  mi_saveas("D:\\FEMM Source\\femmTestFiles\\probe_passB.fem")
  mi_probdef(f, "millimeters", "planar", 1e-008, 35, 30)
  mi_modifymaterial("Amorphous",     5, 0)  -- Cduct (Sigma)
  mi_modifymaterial("Amorphous gap", 5, 0)
  mi_modifymaterial("Amorphous",     6, 0)  -- Lam_d
  mi_modifymaterial("Amorphous gap", 6, 0)
  mi_analyse(1)
  mi_loadsolution()
  arrB_all = {}
  xi = 1
  while xi <= NX do
    arrB_all[xi] = sample_line(X_LIST[xi])
    xi = xi + 1
  end
  mo_close()
  mi_close()

  xi = 1
  while xi <= NX do
    arrA = arrA_all[xi]
    arrB = arrB_all[xi]
    compare(arrA, arrB, format("f=%d Hz, x=%.3f mm", f, X_LIST[xi]))
    -- Dump per-point CSV for plotting (Bx/By are FEMM components)
    local k = 1
    while k <= NPTS do
      write(cf, format("%d,%.6f,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
        f, arrA[k].x, arrA[k].y,
        abs(arrA[k].A),  abs(arrB[k].A),
        abs(arrA[k].B1), abs(arrB[k].B1),
        abs(arrA[k].B2), abs(arrB[k].B2)))
      k = k + 1
    end
    xi = xi + 1
  end
  logln("")
  fi = fi + 1
end

closefile(fp)
closefile(cf)
-- Batch mode: avoid blocking modal dialog to ensure files are flushed/closed.
quit()
