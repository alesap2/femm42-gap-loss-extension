-- gap_battery4_case.lua (Lua 4.0 compatible)
-- Battery 4: Wang gap-scaling of perpendicular flux losses.
--
-- Objective: reproduce Wang et al. (2017) finding that fringing losses scale
-- with gap length.  Instead of total P_side, this battery computes P_x (the
-- fraction of P_side attributable to B_x, perpendicular to laminations) via
-- a 2-D grid sample: P_x = f_Bx * P_side, where
--
--   f_Bx = sum(|Bx_i|^2) / sum(|Bx_i|^2 + |By_i|^2)   over the grid
--
-- Under LT0_OFF both B_x and B_y share the same loss coefficient, so f_Bx
-- is the exact fractional contribution of perpendicular-flux losses.
--
-- Two modes are run for each (gap, f, Bn) tuple:
--   LT0_OFF  -- standard tanh model (baseline, used for P_x decomposition)
--   LT2_ON   -- Bessel correction for B_x (comparison)
--
-- Geometry: lower "Amorphous gap" block  y in [14, 34.3] mm (fixed)
--            upper "Amorphous gap" block  y in [36.3+gap_extra, 56.6+gap_extra]
-- Grid: x in {1,3,5,7,9,11} mm (6 pts)
--       y_lower in {15,17,...,33} mm (10 pts)
--       y_upper in {37,39,...,55}+gap_extra mm (10 pts)
-- Total: 6 * 20 = 120 points per case.

cfg_path   = "D:\\FEMM Source\\femmTestFiles\\gap_battery4_case.cfg"
trace_path = "D:\\FEMM Source\\femmTestFiles\\gap_battery4_out\\gap_battery4_trace.log"

-- ── Utilities ─────────────────────────────────────────────────────────────────

function trace(s)
  local fp = openfile(trace_path, "a")
  if fp ~= nil then
    write(fp, s, "\n")
    closefile(fp)
  end
end

function trim(s)
  if s == nil then return "" end
  s = gsub(s, "^%s+", "")
  s = gsub(s, "%s+$", "")
  return s
end

function read_cfg(path)
  local t = {}
  local fp = openfile(path, "r")
  if fp == nil then return nil end
  local line = read(fp, "*l")
  while line ~= nil do
    local p = strfind(line, "=")
    if p ~= nil then
      local k = trim(strsub(line, 1, p-1))
      local v = trim(strsub(line, p+1))
      t[k] = v
    end
    line = read(fp, "*l")
  end
  closefile(fp)
  return t
end

function read_all(path)
  local fp = openfile(path, "r")
  if fp == nil then return nil end
  local txt = ""
  local line = read(fp, "*l")
  while line ~= nil do
    txt = txt .. line .. "\n"
    line = read(fp, "*l")
  end
  closefile(fp)
  return txt
end

function write_all(path, txt)
  local fp = openfile(path, "w")
  if fp == nil then return 0 end
  write(fp, txt)
  closefile(fp)
  return 1
end

function patch_tag(blk, tag, val)
  local patt = "<" .. tag .. "> = [^\n]+"
  if strfind(blk, "<" .. tag .. ">") ~= nil then
    blk = gsub(blk, patt, "<" .. tag .. "> = " .. val, 1)
  else
    blk = gsub(blk, "\n  <EndBlock>", "\n    <" .. tag .. "> = " .. val .. "\n  <EndBlock>", 1)
  end
  return blk
end

function patch_material_and_save(src_fem, tmp_fem, lam_type, perp_enable, d_lam_str)
  local txt = read_all(src_fem)
  if txt == nil then return 0 end

  -- 1. Patch "Amorphous" block (d_lam only).
  local a1, a2 = strfind(txt, "\"Amorphous\"")
  if a1 ~= nil then
    local b1, b2 = strfind(txt, "<EndBlock>", a2)
    if b1 ~= nil then
      local ablk = strsub(txt, a2, b2)
      ablk = patch_tag(ablk, "d_lam", d_lam_str)
      txt = strsub(txt, 1, a2-1) .. ablk .. strsub(txt, b2+1)
    end
  end

  -- 2. Patch "Amorphous gap" block.
  local i1, i2 = strfind(txt, "\"Amorphous gap\"")
  if i1 == nil then return 0 end
  local j1, j2 = strfind(txt, "<EndBlock>", i2)
  if j1 == nil then return 0 end

  local blk = strsub(txt, i2, j2)
  blk = patch_tag(blk, "LamType",   lam_type)
  blk = patch_tag(blk, "PerpLenz",  perp_enable)
  blk = patch_tag(blk, "d_lam",     d_lam_str)

  if tonumber(perp_enable) ~= 0 then
    local sig_start, sig_end = strfind(blk, "<Sigma> = [^\n]+")
    if sig_start ~= nil then
      local sig_line = strsub(blk, sig_start, sig_end)
      blk = patch_tag(blk, "Sigma_t", strsub(sig_line, 11))
    end
  end

  txt = strsub(txt, 1, i2-1) .. blk .. strsub(txt, j2+1)
  return write_all(tmp_fem, txt)
end

function mean_abs_bn_line(x0, y0, x1, y1)
  mo_clearcontour()
  mo_addcontour(x0, y0)
  mo_addcontour(x1, y1)
  local flux_wb, bn_avg = mo_lineintegral(0)
  if bn_avg == nil then bn_avg = 0 end
  local mx = (x0 + x1) * 0.5
  local my = (y0 + y1) * 0.5
  local Av, Bxv, Byv = mo_getpointvalues(mx, my)
  if Bxv == nil then Bxv = 0 end
  if Byv == nil then Byv = 0 end
  return abs(bn_avg), abs(Bxv), abs(Byv)
end

function append_line(path, line)
  local fp = openfile(path, "a")
  if fp ~= nil then
    write(fp, line, "\n")
    closefile(fp)
  end
end

-- ── Read config ───────────────────────────────────────────────────────────────

cfg = read_cfg(cfg_path)
if cfg == nil then trace("cfg=nil") ; quit() end
trace("cfg ok")

case_id      = cfg["case_id"]
src_fem      = cfg["src_fem"]
tmp_fem      = cfg["tmp_fem"]
summary_csv  = cfg["summary_csv"]
mode_name    = cfg["mode_name"]

freq_hz      = tonumber(cfg["freq_hz"])
gap_mm       = tonumber(cfg["gap_mm"])
base_gap_mm  = tonumber(cfg["base_gap_mm"])
target_bn_t  = tonumber(cfg["target_bn_t"])
seed_i_a     = tonumber(cfg["seed_i_a"])
lam_type     = tonumber(cfg["lam_type"])
perp_enable  = tonumber(cfg["perp_enable"])
d_lam_str    = cfg["d_lam_mm"]

h_x0  = tonumber(cfg["h_x0"]);  h_y0 = tonumber(cfg["h_y0"])
h_x1  = tonumber(cfg["h_x1"]);  h_y1 = tonumber(cfg["h_y1"])

gap_extra = gap_mm - base_gap_mm

-- ── Patch and open model ──────────────────────────────────────────────────────

ok = patch_material_and_save(src_fem, tmp_fem, lam_type, perp_enable, d_lam_str)
if ok == 0 then trace("patch failed") ; quit() end
trace("patch ok d_lam=" .. d_lam_str)

open(tmp_fem)
mi_selectgroup(5)
mi_movetranslate(0, gap_extra)
mi_clearselected()
mi_probdef(freq_hz, "millimeters", "planar", 1e-008, 35, 30)
trace("setup ok")

-- ── Calibration solve (1 A) ───────────────────────────────────────────────────

mi_modifycircprop("coil", 1, 1.0)
mi_analyze(1)
mi_loadsolution()
local cal_ci, cal_vi, cal_fi = mo_getcircuitproperties("coil")
local L_cal = 0
if abs(cal_ci) > 1e-18 then
  L_cal = abs(cal_fi) / abs(cal_ci)
end
mo_close()
trace("cal ok L_cal=" .. L_cal)

-- Tune current for target_bn_t.
local N_turns = 8
local S_m2 = 14e-3 * 35e-3
tuned_i = seed_i_a
if target_bn_t > 0 and L_cal > 1e-20 then
  tuned_i = target_bn_t * N_turns * S_m2 / L_cal
end

-- ── Final solve ───────────────────────────────────────────────────────────────

mi_modifycircprop("coil", 1, tuned_i)
mi_analyze(1)
mi_loadsolution()
trace("final solve ok")

bn_mean, bx_h_mean, by_h_mean = mean_abs_bn_line(h_x0, h_y0, h_x1, h_y1)

-- Total losses.
mo_groupselectblock(0)
P_tot = mo_blockintegral(6)
mo_clearblock()

-- Core region losses (all blocks).
mo_selectblock(38.7, 23.3)
mo_selectblock(37.6, 45.7 + gap_extra)
mo_selectblock(4.4, 17.15)
mo_selectblock(4.4, 53.45 + gap_extra)
mo_selectblock(22.0, 35.3)
mo_selectblock(-8.8, 35.3)
P_core = mo_blockintegral(6)
mo_clearblock()

-- "Amorphous gap" blocks only (inner-leg fringing zone).
mo_selectblock(4.4, 17.15)
mo_selectblock(4.4, 53.45 + gap_extra)
P_side = mo_blockintegral(6)
mo_clearblock()

ci, vi, fi = mo_getcircuitproperties("coil")
Zmag = 0; Lmag = 0
if abs(ci) > 1e-18 then
  Zmag = abs(vi) / abs(ci)
  Lmag = abs(fi) / abs(ci)
end

-- ── f_Bx grid: 2-D field sampling in both "Amorphous gap" blocks ─────────────
-- x grid: 6 points at 1, 3, 5, 7, 9, 11 mm (inside inner leg width [0,14])
-- y lower block [14, 34.3]: 10 points at 15, 17, ..., 33 mm (step 2)
-- y upper block [36.3+gap_extra, 56.6+gap_extra]: 10 pts at 37..55+gap_extra
--
-- Under LT0_OFF: loss_coeff_x = loss_coeff_y (same tanh formula)
--   => P_x = f_Bx * P_side  exactly.
-- Under LT2_ON: loss_coeff_x (Bessel) != loss_coeff_y (tanh)
--   => f_Bx * P_side is only an approximation; used for visual comparison.
-- The primary P_x extraction is done on LT0_OFF in Python post-processing.

local sum_bx2 = 0.0
local sum_b2  = 0.0
local grid_n  = 0

-- x-sample positions (mm)
local gx = {}
gx[1]=1.0; gx[2]=3.0; gx[3]=5.0; gx[4]=7.0; gx[5]=9.0; gx[6]=11.0
local ngx = 6

-- y lower-block sample positions (mm)
local gyl = {}
gyl[1]=15.0; gyl[2]=17.0; gyl[3]=19.0; gyl[4]=21.0; gyl[5]=23.0
gyl[6]=25.0; gyl[7]=27.0; gyl[8]=29.0; gyl[9]=31.0; gyl[10]=33.0
local ngyl = 10

-- y upper-block sample positions, base values (mm); shifted by gap_extra in loop
local gyh = {}
gyh[1]=37.0; gyh[2]=39.0; gyh[3]=41.0; gyh[4]=43.0; gyh[5]=45.0
gyh[6]=47.0; gyh[7]=49.0; gyh[8]=51.0; gyh[9]=53.0; gyh[10]=55.0
local ngyh = 10

local xi = 1
while xi <= ngx do
  local gxv = gx[xi]

  -- lower block
  local yi = 1
  while yi <= ngyl do
    local gyv = gyl[yi]
    local Av, Bxv, Byv = mo_getpointvalues(gxv, gyv)
    if Bxv == nil then Bxv = 0.0 end
    if Byv == nil then Byv = 0.0 end
    local bxm = abs(Bxv)
    local bym = abs(Byv)
    sum_bx2 = sum_bx2 + bxm * bxm
    sum_b2  = sum_b2  + bxm * bxm + bym * bym
    grid_n  = grid_n + 1
    yi = yi + 1
  end

  -- upper block (shifted by gap_extra)
  yi = 1
  while yi <= ngyh do
    local gyv = gyh[yi] + gap_extra
    local Av, Bxv, Byv = mo_getpointvalues(gxv, gyv)
    if Bxv == nil then Bxv = 0.0 end
    if Byv == nil then Byv = 0.0 end
    local bxm = abs(Bxv)
    local bym = abs(Byv)
    sum_bx2 = sum_bx2 + bxm * bxm
    sum_b2  = sum_b2  + bxm * bxm + bym * bym
    grid_n  = grid_n + 1
    yi = yi + 1
  end

  xi = xi + 1
end

local f_bx = 0.0
if sum_b2 > 1e-30 then
  f_bx = sum_bx2 / sum_b2
end
trace("grid ok n=" .. grid_n .. " f_bx=" .. format("%.6f", f_bx))

-- ── Summary row ───────────────────────────────────────────────────────────────
-- Columns: case_id, mode, gap_mm, freq_hz, target_bn_t, d_lam_mm,
--          lam_type, perp_enable, current_a,
--          bn_mean_t, bx_h_mean_t, by_h_mean_t,
--          l_h, z_ohm, p_tot_w, p_core_w, p_side_w,
--          f_bx, grid_n

summary_line = format(
  "%s,%s,%.6f,%d,%.6f,%s,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.8f,%d",
  case_id, mode_name, gap_mm, freq_hz, target_bn_t, d_lam_str,
  lam_type, perp_enable,
  tuned_i, bn_mean, bx_h_mean, by_h_mean,
  Lmag, Zmag, P_tot, P_core, P_side,
  f_bx, grid_n)
append_line(summary_csv, summary_line)
trace("summary ok")

mo_close()
mi_close()
trace("done")
quit()
