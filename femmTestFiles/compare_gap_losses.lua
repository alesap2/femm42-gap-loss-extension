-- compare_gap_losses.lua  (Lua 4.0 / FEMM 4.2)
-- ON vs OFF: total losses & losses integrated over "Amorphous gap" blocks
src     = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev2.fem"
tmp_on  = "D:\\FEMM Source\\femmTestFiles\\_gap_on.fem"
tmp_off = "D:\\FEMM Source\\femmTestFiles\\_gap_off.fem"
outtxt  = "D:\\FEMM Source\\femmTestFiles\\compare_gap_results.txt"

buf = ""
function blog(s) buf = buf .. s .. "\n" end

function run_case(tag, enable_flag, tmppath)
  open(src)
  mi_saveas(tmppath)
  if enable_flag == 1 then
    mi_setmatperplenz("Amorphous gap", 1, 1)
  else
    mi_setmatperplenz("Amorphous gap", 0)
  end
  mi_saveas(tmppath)
  mi_analyze(1)
  mi_loadsolution()

  -- Total ohmic + hysteresis losses across the whole problem
  mo_groupselectblock(0)               -- 0 means "all"
  P_tot = mo_blockintegral(6)          -- 6 = total losses
  mo_clearblock()

  -- Losses inside "Amorphous gap" labels only:
  -- select by clicking inside each "Amorphous gap" label centroid.
  -- The label coordinates for "Amorphous gap" (BlockType==3 in 0-based, 4 in this file)
  -- are at x in {7, 37, 22, -8.8}, y=35.3
  mo_clearblock()
  mo_selectblock(7,   35.3)
  mo_selectblock(37,  35.3)
  mo_selectblock(22,  35.3)
  mo_selectblock(-8.8,35.3)
  P_gap = mo_blockintegral(6)
  mo_clearblock()

  ci, vi, fi = mo_getcircuitproperties("coil")

  blog(format("[%s] P_total       = %.4f W", tag, P_tot))
  blog(format("[%s] P_gap_blocks  = %.4f W", tag, P_gap))
  blog(format("[%s] P_rest        = %.4f W", tag, P_tot - P_gap))
  blog(format("[%s] |I|=%.2f A  |V|=%.4f V  |Z|=%.4e Ohm",
              tag, abs(ci), abs(vi), abs(vi)/abs(ci)))
  blog("")

  mo_close()
  mi_close()
  return P_tot, P_gap
end

blog("=== rev2 PerpLenz gap-loss comparison @ 100 kHz ===")
blog("Material 'Amorphous gap' currently has LamType=1 (lam plane XZ, stacked Y)")
blog("")
Pt_on,  Pg_on  = run_case("ON",  1, tmp_on)
Pt_off, Pg_off = run_case("OFF", 0, tmp_off)

blog(format("Delta P_total      = %.4f W   (%.2f %%)",
            Pt_on-Pt_off, 100*(Pt_on-Pt_off)/Pt_off))
blog(format("Delta P_gap_blocks = %.4f W   (%.2f %%)",
            Pg_on-Pg_off, 100*(Pg_on-Pg_off)/Pg_off))

fh = openfile(outtxt, "w")
write(fh, buf)
closefile(fh)
