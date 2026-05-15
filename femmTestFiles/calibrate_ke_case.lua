-- calibrate_ke_case.lua (Lua 4.0 compatible)
-- Single-case eddy-loss calibration: gap=2mm base, LamType=0, PerpLenz=0.
-- Extracts P_side via blockintegral(3) and block area via blockintegral(5).
-- Config file: calibrate_ke_case.cfg

-- Base model: pourleroi_cc_magnetostatic_rev4.fem (no gap, single "Amorphous gap" block at 7,35.3)
cfg_path   = "D:\\FEMM Source\\femmTestFiles\\calibrate_ke_case.cfg"
trace_path = "D:\\FEMM Source\\femmTestFiles\\calibrate_ke_out\\calibrate_ke_trace.log"

function trace(s)
  local fp = openfile(trace_path, "a")
  if fp ~= nil then write(fp, s, "\n"); closefile(fp) end
end

function trim(s)
  if s == nil then return "" end
  s = gsub(s, "^%s+", ""); s = gsub(s, "%s+$", "")
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
      t[trim(strsub(line,1,p-1))] = trim(strsub(line,p+1))
    end
    line = read(fp, "*l")
  end
  closefile(fp)
  return t
end

function read_all(path)
  local fp = openfile(path, "r")
  if fp == nil then return nil end
  local txt = ""; local line = read(fp, "*l")
  while line ~= nil do txt = txt .. line .. "\n"; line = read(fp, "*l") end
  closefile(fp)
  return txt
end

function write_all(path, txt)
  local fp = openfile(path, "w")
  if fp == nil then return 0 end
  write(fp, txt); closefile(fp); return 1
end

-- Patch LamType, PerpLenz, and optionally d_lam in the "Amorphous gap" block.
-- lam_type and perp_enable are always set; d_lam_mm is applied only when > 0.
function patch_tag(blk, tag, val)
  local patt = "<" .. tag .. "> = [^\n]+"
  if strfind(blk, "<" .. tag .. ">") ~= nil then
    blk = gsub(blk, patt, "<" .. tag .. "> = " .. val, 1)
  else
    blk = gsub(blk, "\n  <EndBlock>", "\n    <" .. tag .. "> = " .. val .. "\n  <EndBlock>", 1)
  end
  return blk
end

function patch_and_save(src, dst, lam_type, perp_enable, d_lam_mm)
  local txt = read_all(src)
  if txt == nil then return 0 end
  local i1, i2 = strfind(txt, "\"Amorphous gap\"")
  if i1 == nil then return 0 end
  local j1, j2 = strfind(txt, "<EndBlock>", i2)
  if j1 == nil then return 0 end
  local blk = strsub(txt, i2, j2)
  blk = patch_tag(blk, "LamType", lam_type)
  blk = patch_tag(blk, "PerpLenz", perp_enable)
  if d_lam_mm ~= nil and d_lam_mm > 0 then
    blk = patch_tag(blk, "d_lam", format("%.6f", d_lam_mm))
  end
  txt = strsub(txt, 1, i2-1) .. blk .. strsub(txt, j2+1)
  return write_all(dst, txt)
end

function mean_abs_bn_line(x0, y0, x1, y1)
  mo_clearcontour()
  mo_addcontour(x0, y0); mo_addcontour(x1, y1)
  local flux_wb, bn_avg = mo_lineintegral(0)
  if bn_avg == nil then bn_avg = 0 end
  return abs(bn_avg)
end

function append_csv(path, line)
  local fp = openfile(path, "a")
  if fp ~= nil then write(fp, line, "\n"); closefile(fp) end
end

-- --------------------------------------------------------------------------
cfg = read_cfg(cfg_path)
if cfg == nil then trace("cfg=nil"); quit() end

case_id     = cfg["case_id"]
src_fem     = cfg["src_fem"]
tmp_fem     = cfg["tmp_fem"]
summary_csv = cfg["summary_csv"]
freq_hz     = tonumber(cfg["freq_hz"])
target_bn   = tonumber(cfg["target_bn_t"])
seed_i      = tonumber(cfg["seed_i_a"])
h_x0        = tonumber(cfg["h_x0"])
h_y0        = tonumber(cfg["h_y0"])
h_x1        = tonumber(cfg["h_x1"])
h_y1        = tonumber(cfg["h_y1"])
-- d_lam in mm; 0 means keep the value already in the .fem file
local d_lam_mm = 0
if cfg["d_lam_mm"] ~= nil then d_lam_mm = tonumber(cfg["d_lam_mm"]) end

ok = patch_and_save(src_fem, tmp_fem, 0, 0, d_lam_mm)
if ok == 0 then trace("patch failed"); quit() end

-- rev4 has no gap — open directly, no geometry move needed.
open(tmp_fem)
mi_probdef(freq_hz, "millimeters", "planar", 1e-008, 35, 30)

-- Calibration pass at 1 A.
mi_modifycircprop("coil", 1, 1.0)
mi_analyze(1); mi_loadsolution()
local cal_ci, cal_vi, cal_fi = mo_getcircuitproperties("coil")
local L_cal = 0
if abs(cal_ci) > 1e-18 then L_cal = abs(cal_fi) / abs(cal_ci) end
mo_close()
trace("L_cal=" .. L_cal)

-- Tuned current: I = B_n * N * S / L
local N = 8; local S = 14e-3 * 35e-3
local tuned_i = seed_i
if target_bn > 0 and L_cal > 1e-20 then
  tuned_i = target_bn * N * S / L_cal
end

-- Final solve.
mi_modifycircprop("coil", 1, tuned_i)
mi_analyze(1); mi_loadsolution()
trace("final solve ok")

-- B_n mean on horizontal calibration line.
local bn_mean = mean_abs_bn_line(h_x0, h_y0, h_x1, h_y1)

-- Single "Amorphous gap" block at (7, 35.3) — rev4, no gap.
mo_selectblock(7, 35.3)
local P_side3 = mo_blockintegral(3)   -- eddy losses [W] via Im(mu_eff)
local A_side  = mo_blockintegral(5)   -- cross-section area [m²]
mo_clearblock()

-- Full core eddy losses (all Amorphous + Amorphous gap blocks).
mo_selectblock(7, 35.3)      -- Amorphous gap (inner leg center)
mo_selectblock(4.4, 17.15)   -- Amorphous (inner leg bottom-left)
mo_selectblock(4.4, 83.45)   -- Amorphous (inner leg top-left)
mo_selectblock(38.7, 23.3)   -- Amorphous (right leg bottom)
mo_selectblock(37, 35.3)     -- Amorphous (right leg center)
mo_selectblock(37.6, 75.7)   -- Amorphous (right leg top)
mo_selectblock(21.9, 3.8)    -- Amorphous (bottom yoke)
mo_selectblock(21.2, 96.3)   -- Amorphous (top yoke)
local P_core3 = mo_blockintegral(3)
mo_clearblock()

local line = format("%s,%d,%.6f,%.6e,%.6e,%.6e,%.6e",
  case_id, freq_hz, target_bn, tuned_i, bn_mean, P_side3, A_side)
-- Also append P_core3 as extra column.
line = line .. format(",%.6e", P_core3)
append_csv(summary_csv, line)
trace("done: " .. line)

mo_close(); mi_close(); quit()
