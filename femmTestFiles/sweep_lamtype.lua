-- sweep_lamtype.lua  (Lua 4.0 / FEMM 4.2)
-- Test "Amorphous gap" with all 3 LamType options x 3 PerpLenz cases.
-- Reports: total losses, losses on side gaps, flux linkage, |Z|.
src     = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev2.fem"
legacy  = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev2_legacy.fem"
tmpfem  = "D:\\FEMM Source\\femmTestFiles\\_sweep_tmp.fem"
outtxt  = "D:\\FEMM Source\\femmTestFiles\\sweep_lamtype_results.txt"

buf = ""
function blog(s) buf = buf .. s .. "\n" end

-- Run analysis on a .fem already saved
function analyze(tag)
  mi_analyze(1)
  mi_loadsolution()

  mo_groupselectblock(0)
  P_tot = mo_blockintegral(6)
  mo_clearblock()

  -- CORE-ONLY losses: select all non-copper, non-air blocks
  -- "Amorphous" labels (BlockType=1):
  mo_selectblock(38.7, 23.3)
  mo_selectblock(37.6, 50.7)
  -- "Amorphous gap" labels (BlockType=3):
  mo_selectblock(4.4, 17.15)
  mo_selectblock(4.4, 58.45)
  mo_selectblock(22.0, 35.3)
  mo_selectblock(-8.8, 35.3)
  -- u1..u7 labels (BlockType=5..11):
  mo_selectblock(231.2, 76.9)
  mo_selectblock(221.9, 118.1)
  mo_selectblock(204.4, 157.2)
  mo_selectblock(179.3, 192.6)
  mo_selectblock(147.2, 222.7)
  mo_selectblock(109.4, 246.4)
  mo_selectblock(67.2, 262.3)
  P_core = mo_blockintegral(6)
  mo_clearblock()

  -- side gap labels (Amorphous gap columns, x=4.4)
  mo_selectblock(4.4, 17.15)
  mo_selectblock(4.4, 58.45)
  P_side = mo_blockintegral(6)
  mo_clearblock()

  ci, vi, fi = mo_getcircuitproperties("coil")

  A, B1, B2 = mo_getpointvalues(4.4, 17.15)
  Bx_in = abs(B1); By_in = abs(B2)

  blog(format("[%s] P_tot=%.3f W  P_core=%.3f W  P_side=%.3f W  |flx|=%.2e Wb  |Z|=%.3e",
              tag, P_tot, P_core, P_side, abs(fi), abs(vi)/abs(ci)))
  blog(format("        B inside side gap label (4.4,17.15): |Bx|=%.3e  |By|=%.3e T",
              Bx_in, By_in))
  mo_close()
  mi_close()
  return P_tot, P_side, abs(fi)
end

-- Set LamType + PerpLenz for "Amorphous gap" by patching the .fem text.
-- (mi_modifymaterial does not expose LamType in the public API consistently,
--  so we directly edit the saved .fem then reopen.)
function patch_fem(lam_type, perp_enable, perp_model)
  -- read source
  fp = openfile(src, "r")
  txt = ""
  line = read(fp, "*l")
  while line ~= nil do
    txt = txt .. line .. "\n"
    line = read(fp, "*l")
  end
  closefile(fp)

  -- locate "Amorphous gap" block and patch LamType
  -- naive single replacement: only first occurrence of LamType inside that block
  i1, i2 = strfind(txt, "\"Amorphous gap\"")
  i3, i4 = strfind(txt, "<EndBlock>", i2)
  -- inside [i2, i3] is the Amorphous gap block content
  blk = strsub(txt, i2, i3)
  -- replace LamType = X
  blk = gsub(blk, "<LamType> = %d+", "<LamType> = " .. lam_type, 1)
  -- replace or insert PerpLenz
  if strfind(blk, "<PerpLenz>") then
    blk = gsub(blk, "<PerpLenz> = %d+", "<PerpLenz> = " .. perp_enable, 1)
  else
    blk = gsub(blk, "<sigma_n>",
               "<PerpLenz> = " .. perp_enable .. "\n    <sigma_n>", 1)
  end
  txt = strsub(txt, 1, i2-1) .. blk .. strsub(txt, i3+1)

  -- write tmp
  fo = openfile(tmpfem, "w")
  write(fo, txt)
  closefile(fo)
end

blog("=== Sweep LamType x PerpLenz on 'Amorphous gap' (100 kHz) ===")
blog("")

-- Legacy reference once
open(legacy)
Pleg, Sleg, Fleg = analyze("LEGACY-public")
blog("")

-- Test matrix
cases = {
  {0, 0, 0, "LamType=0 PerpLenz=OFF (isotropic)"},
  {1, 0, 0, "LamType=1 PerpLenz=OFF (lam plane XZ)"},
  {1, 1, 1, "LamType=1 PerpLenz=ON  model=1 (Wperp=d_lam)"},
  {1, 1, 2, "LamType=1 PerpLenz=ON  model=2 (Wperp=bbox)"},
  {2, 0, 0, "LamType=2 PerpLenz=OFF (lam plane YZ)"},
  {2, 1, 1, "LamType=2 PerpLenz=ON  model=1 (Wperp=d_lam)"},
  {2, 1, 2, "LamType=2 PerpLenz=ON  model=2 (Wperp=bbox)"},
}

n = 1
while cases[n] ~= nil do
  c = cases[n]
  patch_fem(c[1], c[2], c[3])
  open(tmpfem)
  analyze(c[4])
  blog("")
  n = n + 1
end

fh = openfile(outtxt, "w")
write(fh, buf)
closefile(fh)
