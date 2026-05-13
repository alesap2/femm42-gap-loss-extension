-- probe_perp_lenz.lua  (Lua 4.0 compatible)
-- Validation VP-1, VP-2, VP-5 for perpendicular Lenz mu_perp(omega) model.
--
-- Invoked once per frequency by run_perp_lenz.ps1.
-- Each run does 3 solves (REF, VP-5, VP-2/1) and appends to shared CSV/log.
--
-- Material: "Amorphous gap"  (LamType==2, sigma_t>0)  => perp slot = Bx
-- Geometry: pourleroi_cc_magnetostatic_rev3.fem (planar, mm)
-- Probe: x=7mm, y in [14.0, 60.5] mm

src  = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev3.fem"
logp = "D:\\FEMM Source\\femmTestFiles\\probe_perp_lenz_results.txt"
csvp = "D:\\FEMM Source\\femmTestFiles\\probe_perp_lenz_data.csv"
idxp = "D:\\FEMM Source\\femmTestFiles\\probe_perp_lenz_idx.txt"
tmp  = "D:\\FEMM Source\\femmTestFiles\\probe_pl_tmp.fem"
MAT  = "Amorphous gap"

NPTS   = 201
Y0     = 14.0
Y1     = 60.5
DY     = (Y1 - Y0) / (NPTS - 1)
XPROBE = 7.0
FREQS  = {1000, 10000, 100000}

-- Read frequency index.
fi = 1
h = openfile(idxp, "r")
if h ~= nil then
  sv = read(h, "*l")
  closefile(h)
  if sv ~= nil then fi = tonumber(sv) end
end
h2 = openfile(idxp, "w")
write(h2, tostring(fi + 1), "\n")
closefile(h2)

f = FREQS[fi]
if f == nil then quit() end

first_run = (fi == 1)

-- CSV data is written inside sample_line().
-- Log text is accumulated in logbuf and written AFTER all femm work.
logbuf = ""

function logadd(s)
  logbuf = logbuf .. s .. "\n"
end

-- Accumulate CSV rows in csvbuf; written after all femm work.
csvbuf = ""

function sample_line(tag, freq)
  mbx = 0
  mby = 0
  ii = 1
  while ii <= NPTS do
    yy = Y0 + (ii - 1) * DY
    A, B1, B2 = mo_getpointvalues(XPROBE, yy)
    if B1 ~= nil and B2 ~= nil then
      bx = abs(B1)
      by = abs(B2)
      csvbuf = csvbuf .. tag .. "," .. tostring(freq) .. ","
                      .. format("%.4f", yy) .. ","
                      .. format("%.6e", bx) .. ","
                      .. format("%.6e", by) .. "\n"
      if bx > mbx then mbx = bx end
      if by > mby then mby = by end
    end
    ii = ii + 1
  end
  return mbx, mby
end

logadd(format("---- f = %d Hz ----", f))

-- REF
open(src)
mi_saveas(tmp)
mi_probdef(f, "millimeters", "planar", 1e-008, 35, 30)
mi_analyse(1)
mi_loadsolution()
bxR, byR = sample_line("ref", f)
mo_close()
mi_close()
logadd(format("  REF    max|Bx|=%.4e T  max|By|=%.4e T", bxR, byR))

-- VP-5: bPerpLenz disabled
open(src)
mi_saveas(tmp)
mi_probdef(f, "millimeters", "planar", 1e-008, 35, 30)
mi_setmatperplenz(MAT, 0)
mi_analyse(1)
mi_loadsolution()
bxV5, byV5 = sample_line("vp5", f)
mo_close()
mi_close()
p5 = 0
if bxR > 0 then p5 = 100 * abs(bxV5 - bxR) / bxR end
logadd(format("  VP-5   max|Bx|=%.4e T  max|By|=%.4e T  diff_Bx=%.5f%%", bxV5, byV5, p5))

-- VP-2/VP-1: bPerpLenz enabled
open(src)
mi_saveas(tmp)
mi_probdef(f, "millimeters", "planar", 1e-008, 35, 30)
mi_setmatperplenz(MAT, 1)
mi_analyse(1)
mi_loadsolution()
bxV2, byV2 = sample_line("vp2", f)
mo_close()
mi_close()
pbx = 0
pby = 0
if bxR > 0 then pbx = 100 * abs(bxV2 - bxR) / bxR end
if byR > 0 then pby = 100 * abs(byV2 - byR) / byR end
logadd(format("  VP-2/1 max|Bx|=%.4e T  max|By|=%.4e T  dBx=%.2f%%  dBy=%.2f%%",
              bxV2, byV2, pbx, pby))
logadd("")

-- Write CSV (all femm work done, file handles now safe).
if first_run then
  cf = openfile(csvp, "w")
  write(cf, "test,f_Hz,y_mm,absBx,absBy\n")
else
  cf = openfile(csvp, "a")
end
write(cf, csvbuf)
closefile(cf)

-- Write log.
if first_run then
  fp = openfile(logp, "w")
  write(fp, "=== probe_perp_lenz.lua ===\n")
  write(fp, "Material: ", MAT, "\n")
  write(fp, format("Probe x=%.1f mm, y=[%.1f,%.1f] mm, %d pts\n\n", XPROBE, Y0, Y1, NPTS))
else
  fp = openfile(logp, "a")
end
write(fp, logbuf)
closefile(fp)

quit()
