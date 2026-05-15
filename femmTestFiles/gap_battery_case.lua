-- gap_battery_case.lua (Lua 4.0 compatible)
-- Executes one FEMM case defined in gap_battery_case.cfg
-- Base .fem is never modified; a patched temporary model is created each run.

cfg_path = "D:\\FEMM Source\\femmTestFiles\\gap_battery_case.cfg"
trace_path = "D:\\FEMM Source\\femmTestFiles\\gap_battery_out\\gap_battery_trace.log"

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
  if fp == nil then
    return nil
  end
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

function patch_material_and_save(src_fem, tmp_fem, lam_type, perp_enable)
  local txt = read_all(src_fem)
  if txt == nil then return 0 end

  local i1, i2 = strfind(txt, "\"Amorphous gap\"")
  if i1 == nil then return 0 end

  local j1, j2 = strfind(txt, "<EndBlock>", i2)
  if j1 == nil then return 0 end

  local blk = strsub(txt, i2, j2)
  blk = patch_tag(blk, "LamType", lam_type)
  blk = patch_tag(blk, "PerpLenz", perp_enable)
  
  -- When PerpLenz is enabled, also set Cduct_t to the same value as Sigma
  -- so the Bessel perpendicular model can compute frequency-dependent eddy losses.
  if tonumber(perp_enable) ~= 0 then
    -- Extract Sigma value from block
    local sig_start, sig_end = strfind(blk, "<Sigma> = [^\n]+")
    if sig_start ~= nil then
      local sig_line = strsub(blk, sig_start, sig_end)
      -- Patch or add Sigma_t to match Sigma (for use as Cduct_t in solver/postproc)
      blk = patch_tag(blk, "Sigma_t", strsub(sig_line, 11))  -- Extract the value part
    end
  end

  txt = strsub(txt, 1, i2-1) .. blk .. strsub(txt, j2+1)
  return write_all(tmp_fem, txt)
end

function mean_abs_bn_line(x0, y0, x1, y1, npts)
  -- FEMM native B.n line integral (exact, uses FEM shape functions).
  -- mo_lineintegral(0) returns TWO values:
  --   [1] flux_wb  = (A_start - A_end)*Depth  [Wb, complex phasor]
  --   [2] bn_avg   = flux_wb / (length_m * Depth)  [T, complex phasor]
  -- abs(bn_avg) = peak amplitude of average B.n.
  mo_clearcontour()
  mo_addcontour(x0, y0)
  mo_addcontour(x1, y1)
  local flux_wb, bn_avg = mo_lineintegral(0)
  if bn_avg == nil then bn_avg = 0 end
  local bn_mean = abs(bn_avg)
  -- Midpoint sample for Bx/By diagnostics only.
  local mx = (x0 + x1) * 0.5
  local my = (y0 + y1) * 0.5
  local A, Bx, By = mo_getpointvalues(mx, my)
  if Bx == nil then Bx = 0 end
  if By == nil then By = 0 end
  return bn_mean, abs(Bx), abs(By)
end

function append_summary(path, line)
  local fp = openfile(path, "a")
  if fp ~= nil then
    write(fp, line, "\n")
    closefile(fp)
  end
end

function append_profile(path, line)
  local fp = openfile(path, "a")
  if fp ~= nil then
    write(fp, line, "\n")
    closefile(fp)
  end
end

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
mode_name     = cfg["mode_name"]

freq_hz       = tonumber(cfg["freq_hz"])
gap_mm        = tonumber(cfg["gap_mm"])
base_gap_mm   = tonumber(cfg["base_gap_mm"])
target_bn_t   = tonumber(cfg["target_bn_t"])
seed_i_a      = tonumber(cfg["seed_i_a"])
lam_type      = tonumber(cfg["lam_type"])
perp_enable   = tonumber(cfg["perp_enable"])

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

ok = patch_material_and_save(src_fem, tmp_fem, lam_type, perp_enable)
if ok == 0 then
  trace("patch failed")
  quit()
end

trace("patch ok")

open(tmp_fem)
trace("open tmpfem ok")
mi_selectgroup(5)
mi_movetranslate(0, gap_extra)
mi_clearselected()

mi_probdef(freq_hz, "millimeters", "planar", 1e-008, 35, 30)
trace("setup model ok")

-- Calibration pass at 1 A to extract inductance L.
-- Then: I_target = target_Bn * N_turns * S_m2 / L
-- N=8 turns, S=14mm x 35mm (problem depth) = 4.9e-4 m^2.
mi_modifycircprop("coil", 1, 1.0)
mi_analyze(1)
mi_loadsolution()
local cal_ci, cal_vi, cal_fi = mo_getcircuitproperties("coil")
local L_cal = 0
if abs(cal_ci) > 1e-18 then
  L_cal = abs(cal_fi) / abs(cal_ci)
end
mo_close()
trace("cal solve ok, L_cal=" .. L_cal)

local N_turns = 8
local S_m2 = 14e-3 * 35e-3
tuned_i = seed_i_a
if target_bn_t > 0 and L_cal > 1e-20 then
  tuned_i = target_bn_t * N_turns * S_m2 / L_cal
end

-- Final pass (tuned current).
mi_modifycircprop("coil", 1, tuned_i)
mi_analyze(1)
mi_loadsolution()
trace("final solve ok")

bn_mean, bx_h_mean, by_h_mean = mean_abs_bn_line(h_x0, h_y0, h_x1, h_y1, h_n)

mo_groupselectblock(0)
P_tot = mo_blockintegral(6)
mo_clearblock()

-- Core region losses (Amorphous + Amorphous gap labels only).
mo_selectblock(38.7, 23.3)
mo_selectblock(37.6, 45.7 + gap_extra)
mo_selectblock(4.4, 17.15)
mo_selectblock(4.4, 53.45 + gap_extra)
mo_selectblock(22.0, 35.3)
mo_selectblock(-8.8, 35.3)
P_core = mo_blockintegral(6)
mo_clearblock()

-- Side-gap only losses.
mo_selectblock(4.4, 17.15)
mo_selectblock(4.4, 53.45 + gap_extra)
P_side = mo_blockintegral(6)
mo_clearblock()

-- Thin-lam analytical integral over gap pieces.
P_thin31 = 0

ci, vi, fi = mo_getcircuitproperties("coil")
Zmag = 0
Lmag = 0
if abs(ci) > 1e-18 then
  Zmag = abs(vi) / abs(ci)
  Lmag = abs(fi) / abs(ci)
end

summary_line = format("%s,%s,%.6f,%d,%.6f,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e",
  case_id, mode_name, gap_mm, freq_hz, target_bn_t, lam_type, perp_enable,
  tuned_i, bn_mean, bx_h_mean, by_h_mean, Lmag, Zmag, P_tot, P_core, P_side, P_thin31)
append_summary(summary_csv, summary_line)
trace("summary appended")

-- Vertical Bx profile on requested line.
local i = 0
while i < v_n do
  local t = 0
  if v_n > 1 then t = i / (v_n - 1) end
  local y = v_y0 + (v_y1 - v_y0) * t
  local A, Bx, By = mo_getpointvalues(v_x, y)
  if Bx == nil then Bx = 0 end
  if By == nil then By = 0 end
  local pl = format("%s,%s,%.6f,%d,%.6f,%.6f,%.6e,%.6e",
    case_id, mode_name, gap_mm, freq_hz, target_bn_t, y, abs(Bx), abs(By))
  append_profile(profile_csv, pl)
  i = i + 1
end

mo_close()
mi_close()
trace("done")
quit()
