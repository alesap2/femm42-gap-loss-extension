-- gap_battery2_case.lua (Lua 4.0 compatible)
-- New battery: 2 modes (LT0_OFF, LT2_ON), d_lam sweep, optional 2D grid export.
-- Key difference vs gap_battery_case.lua:
--   * Patches BOTH "Amorphous" and "Amorphous gap" blocks with d_lam from cfg.
--   * If grid_enabled=1, exports |Bx|/|By| on a 14x20 grid over the lower P_side block.
-- Base .fem is never modified; a patched temporary model is created each run.

cfg_path = "D:\\FEMM Source\\femmTestFiles\\gap_battery2_case.cfg"
trace_path = "D:\\FEMM Source\\femmTestFiles\\gap_battery2_out\\gap_battery2_trace.log"

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

-- Patch both "Amorphous" and "Amorphous gap" blocks.
-- "Amorphous" block: only d_lam updated (LamType stays 0, LamFill stays as-is).
-- "Amorphous gap" block: LamType, PerpLenz, Sigma_t, d_lam all updated.
-- Order: patch "Amorphous" first (appears earlier in file) to avoid position shift issues
-- when later searching for "Amorphous gap".
function patch_material_and_save(src_fem, tmp_fem, lam_type, perp_enable, d_lam_str)
  local txt = read_all(src_fem)
  if txt == nil then return 0 end

  -- 1. Patch "Amorphous" block (d_lam only).
  --    Search for the exact string "Amorphous" with surrounding double-quotes.
  --    This does NOT match "Amorphous gap" because that string contains ' gap' before
  --    the closing quote, so the substring "Amorphous" (with closing quote) is absent.
  local a1, a2 = strfind(txt, "\"Amorphous\"")
  if a1 ~= nil then
    local b1, b2 = strfind(txt, "<EndBlock>", a2)
    if b1 ~= nil then
      local ablk = strsub(txt, a2, b2)
      ablk = patch_tag(ablk, "d_lam", d_lam_str)
      txt = strsub(txt, 1, a2-1) .. ablk .. strsub(txt, b2+1)
    end
  end

  -- 2. Patch "Amorphous gap" block (re-search after string may have changed length).
  local i1, i2 = strfind(txt, "\"Amorphous gap\"")
  if i1 == nil then return 0 end
  local j1, j2 = strfind(txt, "<EndBlock>", i2)
  if j1 == nil then return 0 end

  local blk = strsub(txt, i2, j2)
  blk = patch_tag(blk, "LamType", lam_type)
  blk = patch_tag(blk, "PerpLenz", perp_enable)
  blk = patch_tag(blk, "d_lam", d_lam_str)

  -- When PerpLenz is enabled, also set Sigma_t to match Sigma so the Bessel
  -- perpendicular model uses the correct conductivity.
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

function mean_abs_bn_line(x0, y0, x1, y1, npts)
  mo_clearcontour()
  mo_addcontour(x0, y0)
  mo_addcontour(x1, y1)
  local flux_wb, bn_avg = mo_lineintegral(0)
  if bn_avg == nil then bn_avg = 0 end
  local bn_mean = abs(bn_avg)
  local mx = (x0 + x1) * 0.5
  local my = (y0 + y1) * 0.5
  local A, Bx, By = mo_getpointvalues(mx, my)
  if Bx == nil then Bx = 0 end
  if By == nil then By = 0 end
  return bn_mean, abs(Bx), abs(By)
end

function append_line(path, line)
  local fp = openfile(path, "a")
  if fp ~= nil then
    write(fp, line, "\n")
    closefile(fp)
  end
end

-- ── Read config ──────────────────────────────────────────────────────────────

cfg = read_cfg(cfg_path)
if cfg == nil then
  trace("cfg=nil")
  quit()
end
trace("cfg ok")

case_id       = cfg["case_id"]
src_fem       = cfg["src_fem"]
tmp_fem       = cfg["tmp_fem"]
summary_csv   = cfg["summary_csv"]
profile_csv   = cfg["profile_csv"]
grid_csv      = cfg["grid_csv"]
mode_name     = cfg["mode_name"]

freq_hz       = tonumber(cfg["freq_hz"])
gap_mm        = tonumber(cfg["gap_mm"])
base_gap_mm   = tonumber(cfg["base_gap_mm"])
target_bn_t   = tonumber(cfg["target_bn_t"])
seed_i_a      = tonumber(cfg["seed_i_a"])
lam_type      = tonumber(cfg["lam_type"])
perp_enable   = tonumber(cfg["perp_enable"])
d_lam_mm_val  = tonumber(cfg["d_lam_mm"])
d_lam_str     = cfg["d_lam_mm"]   -- keep as string for patch_tag
grid_enabled  = tonumber(cfg["grid_enabled"])
if grid_enabled == nil then grid_enabled = 0 end

h_x0          = tonumber(cfg["h_x0"])
h_y0          = tonumber(cfg["h_y0"])
h_x1          = tonumber(cfg["h_x1"])
h_y1          = tonumber(cfg["h_y1"])
h_n           = tonumber(cfg["h_n"])

v_x           = tonumber(cfg["v_x"])
v_y0          = tonumber(cfg["v_y0"])
v_y1_base     = tonumber(cfg["v_y1_base"])
v_n           = tonumber(cfg["v_n"])

gap_extra = gap_mm - base_gap_mm
v_y1 = v_y1_base + gap_extra

-- ── Patch and open model ─────────────────────────────────────────────────────

ok = patch_material_and_save(src_fem, tmp_fem, lam_type, perp_enable, d_lam_str)
if ok == 0 then
  trace("patch failed")
  quit()
end
trace("patch ok, d_lam=" .. d_lam_str)

open(tmp_fem)
mi_selectgroup(5)
mi_movetranslate(0, gap_extra)
mi_clearselected()
mi_probdef(freq_hz, "millimeters", "planar", 1e-008, 35, 30)
trace("setup ok")

-- ── Calibration solve (1 A) → inductance ─────────────────────────────────────

mi_modifycircprop("coil", 1, 1.0)
mi_analyze(1)
mi_loadsolution()
local cal_ci, cal_vi, cal_fi = mo_getcircuitproperties("coil")
local L_cal = 0
if abs(cal_ci) > 1e-18 then
  L_cal = abs(cal_fi) / abs(cal_ci)
end
mo_close()
trace("cal ok, L_cal=" .. L_cal)

local N_turns = 8
local S_m2 = 14e-3 * 35e-3
tuned_i = seed_i_a
if target_bn_t > 0 and L_cal > 1e-20 then
  tuned_i = target_bn_t * N_turns * S_m2 / L_cal
end

-- ── Final solve (tuned current) ───────────────────────────────────────────────

mi_modifycircprop("coil", 1, tuned_i)
mi_analyze(1)
mi_loadsolution()
trace("final solve ok")

bn_mean, bx_h_mean, by_h_mean = mean_abs_bn_line(h_x0, h_y0, h_x1, h_y1, h_n)

-- Total losses (all blocks with sigma>0 respond to harmonic excitation).
mo_groupselectblock(0)
P_tot = mo_blockintegral(6)
mo_clearblock()

-- Core region losses (Amorphous body + gap blocks + air gap).
-- Air blocks (sigma=0) contribute zero to blockintegral(6).
mo_selectblock(38.7, 23.3)
mo_selectblock(37.6, 45.7 + gap_extra)
mo_selectblock(4.4, 17.15)
mo_selectblock(4.4, 53.45 + gap_extra)
mo_selectblock(22.0, 35.3)
mo_selectblock(-8.8, 35.3)
P_core = mo_blockintegral(6)
mo_clearblock()

-- Side-gap blocks only (inner leg, adjacent to gap faces).
mo_selectblock(4.4, 17.15)
mo_selectblock(4.4, 53.45 + gap_extra)
P_side = mo_blockintegral(6)
mo_clearblock()

P_thin31 = 0

ci, vi, fi = mo_getcircuitproperties("coil")
Zmag = 0
Lmag = 0
if abs(ci) > 1e-18 then
  Zmag = abs(vi) / abs(ci)
  Lmag = abs(fi) / abs(ci)
end

-- ── Write summary row ─────────────────────────────────────────────────────────

summary_line = format(
  "%s,%s,%.6f,%d,%.6f,%s,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e",
  case_id, mode_name, gap_mm, freq_hz, target_bn_t, d_lam_str, lam_type, perp_enable,
  tuned_i, bn_mean, bx_h_mean, by_h_mean, Lmag, Zmag, P_tot, P_core, P_side, P_thin31)
append_line(summary_csv, summary_line)
trace("summary ok")

-- ── Vertical Bx/By profile ────────────────────────────────────────────────────

local i = 0
while i < v_n do
  local t = 0
  if v_n > 1 then t = i / (v_n - 1) end
  local y = v_y0 + (v_y1 - v_y0) * t
  local A, Bx, By = mo_getpointvalues(v_x, y)
  if Bx == nil then Bx = 0 end
  if By == nil then By = 0 end
  local pl = format("%s,%s,%.6f,%d,%.6f,%s,%.6f,%.6e,%.6e",
    case_id, mode_name, gap_mm, freq_hz, target_bn_t, d_lam_str, y, abs(Bx), abs(By))
  append_line(profile_csv, pl)
  i = i + 1
end
trace("profile ok")

-- ── 2D grid export (lower P_side block, fixed position) ───────────────────────
-- Samples the inner leg x=[0.5,13.5]mm × y=[14.5,33.5]mm (14×20=280 pts).
-- Only executed when grid_enabled=1.  The lower P_side block label is at
-- (4.4, 17.15) and covers x≈[0,14]mm, y≈[14,34.3]mm (below the gap face).
-- mo_getpointvalues returns complex phasors; abs() gives peak magnitude.

if grid_enabled == 1 then
  local gx0, gx1, gnx = 0.5, 13.5, 14
  local gy0, gy1, gny = 14.5, 33.5, 20
  local ix = 0
  while ix < gnx do
    local gx = gx0
    if gnx > 1 then gx = gx0 + (gx1-gx0) * ix / (gnx-1) end
    local iy = 0
    while iy < gny do
      local gy = gy0
      if gny > 1 then gy = gy0 + (gy1-gy0) * iy / (gny-1) end
      local gA, gBx, gBy = mo_getpointvalues(gx, gy)
      if gBx == nil then gBx = 0 end
      if gBy == nil then gBy = 0 end
      local gl = format("%s,%s,%.6f,%d,%.6f,%s,%.4f,%.4f,%.6e,%.6e",
        case_id, mode_name, gap_mm, freq_hz, target_bn_t, d_lam_str,
        gx, gy, abs(gBx), abs(gBy))
      append_line(grid_csv, gl)
      iy = iy + 1
    end
    ix = ix + 1
  end
  trace("grid ok (280 pts)")
end

mo_close()
mi_close()
trace("done")
quit()
