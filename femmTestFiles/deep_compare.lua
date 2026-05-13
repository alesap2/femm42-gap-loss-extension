-- deep_compare.lua  (Lua 4.0 / FEMM 4.2)
-- Four cases to isolate where the 20% loss difference comes from:
--   A: our solver, PerpLenz=OFF  (should match legacy)
--   B: our solver, PerpLenz=ON   (what user saw: 6.06 W)
--   C: legacy solver (reference: 7.54 W)
-- We also sample mu_eff in the gap by reading B and H at a point.
src     = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev2.fem"
legacy  = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev2_legacy.fem"
tmp_on  = "D:\\FEMM Source\\femmTestFiles\\_dc_on.fem"
tmp_off = "D:\\FEMM Source\\femmTestFiles\\_dc_off.fem"
outtxt  = "D:\\FEMM Source\\femmTestFiles\\deep_compare_results.txt"

buf = ""
function blog(s) buf = buf .. s .. "\n" end

-- gap label centroids in the model
-- "Amorphous gap" block labels (from .fem BlockLabels section, BlockType index 3=Amorphous gap):
-- rows 7,8 in BlockLabels: (4.4,17.15), (4.4,58.45) → Amorphous (LamType=0)
-- rows 9,10: (7,35.3), (37,35.3) → Air
-- rows 11-end: coil turns
-- Actually "Amorphous gap" (index 2, 0-based) labels are the ones with BlockType=2 in file
-- Looking at file: BlockType column is the 3rd field in NumBlockLabels section
-- BlockType=3 → "Amorphous gap" (0-indexed: Air=1, Amorphous=0, AmorGap=2, Litz=3)
-- Wait, block list order: Amorphous=0, Air=1, AmorGap=2, Litz=3, u1=4..u7=10
-- So BlockType=2 = "Amorphous gap"
-- Looking at NumBlockLabels: labels with BlockType=2:
--   7: 4.4, 17.15 (BlockType=3 in file → 0-based=2 = Amorphous gap? No...)
-- Actually in the file, BlockType field = material index 1-based:
-- 1=Amorphous, 2=Air, 3=Amorphous gap, 4=LitzCu, 5=u1...
-- Labels with BlockType=3: (4.4,17.15) and (4.4,58.45)
-- Labels with BlockType=2: (7,35.3) and (37,35.3) → those are Air
-- Wait: let me re-read the NumBlockLabels:
--  231... blocktype=5 (u1)
--  221... blocktype=6 (u2)
--  ...
--  4.4, 17.15 blocktype=3 (Amorphous gap)
--  4.4, 58.45 blocktype=3 (Amorphous gap)  
--  7, 35.3    blocktype=2 (Air)
--  37, 35.3   blocktype=2 (Air)
--  coil turns blocktype=4 (Litz)
-- So the Amorphous gap labels are at x=4.4, y=17.15 and y=58.45
-- And the inner air labels at x=7 and x=37, y=35.3
-- But wait, looking at the image there are horizontal gap blocks in the middle...
-- Let me just probe specific points inside the visible "Amorphous gap" region

-- From the .fem geometry: the Amorphous gap blocks are at y~34.3 (gap region)
-- x=4.4, y=17.15 → inside lower Amorphous gap block
-- x=4.4, y=58.45 → inside upper Amorphous gap block
-- Also x=38.7, y=23.3 and x=37.6, y=50.7 have BlockType=3 (Amorphous)

-- Probe points inside "Amorphous gap" region (horizontal strips at y=34.3±3.5)
-- From geometry: gap is between y=34.3 and y=41.3, width x=0..44
-- Probe at center of gap strip: x=22, y=37.8
XGAP = 22.0
YGAP = 37.8   -- center of horizontal gap strip

function analyze_case(tag, fempath)
  open(fempath)
  mi_analyze(1)
  mi_loadsolution()

  -- Total losses
  mo_groupselectblock(0)
  P_tot = mo_blockintegral(6)
  mo_clearblock()

  -- Losses in blocks at x=4.4 (the two Amorphous gap side columns)
  mo_selectblock(4.4, 17.15)
  mo_selectblock(4.4, 58.45)
  P_side = mo_blockintegral(6)
  mo_clearblock()

  -- Flux linkage and impedance
  ci, vi, fi = mo_getcircuitproperties("coil")

  -- Field at probe point in the gap
  A, B1, B2 = mo_getpointvalues(XGAP, YGAP)
  Bx_gap = abs(B1)
  By_gap = abs(B2)

  -- Also probe inside the solid Amorphous core (not gap)
  A2, Bc1, Bc2 = mo_getpointvalues(4.4, 17.15)
  Bx_core = abs(Bc1)
  By_core = abs(Bc2)

  blog(format("[%s] P_total        = %.4f W", tag, P_tot))
  blog(format("[%s] P_side_gaps     = %.4f W  (x=4.4 col)", tag, P_side))
  blog(format("[%s] |flux linkage|  = %.4e Wb-t", tag, abs(fi)))
  blog(format("[%s] |Z_coil|        = %.4e Ohm", tag, abs(vi)/abs(ci)))
  blog(format("[%s] B at gap center (%.0f,%.1f): |Bx|=%.4e T  |By|=%.4e T",
              tag, XGAP, YGAP, Bx_gap, By_gap))
  blog(format("[%s] B at core (4.4,17.15): |Bx|=%.4e T  |By|=%.4e T",
              tag, Bx_core, By_core))
  blog("")

  mo_close()
  mi_close()
  return P_tot, P_side, abs(fi)
end

function analyze_our(tag, enable_flag, tmppath)
  open(src)
  mi_saveas(tmppath)
  if enable_flag == 1 then
    mi_setmatperplenz("Amorphous gap", 1, 1)
  else
    mi_setmatperplenz("Amorphous gap", 0)
  end
  mi_saveas(tmppath)
  return analyze_case(tag, tmppath)
end

blog("=== Deep comparison: our solver ON vs OFF vs legacy ===")
blog(format("Frequency: 100 kHz   LamType=1 (lam plane XZ, B_y perpendicular)"))
blog(format("d_lam=0.023 mm, sigma_t=0.615 MS/m, LamFill=0.8, BH-curve (mu_init~36600)"))
blog(format("|gamma*a| at 100kHz = d_lam/2 * sqrt(w*mu_init*mu0*sigma_t)"))
blog(format("  = 11.5e-6 * sqrt(6.28e5 * 36600 * 4pi1e-7 * 0.615e6)"))
blog(format("  = 11.5e-6 * ~1.33e5 = ~1.53"))
blog(format("Bessel shape(1.53) ~ 0.75 -> mu_perp ~ 0.8*36600*0.75+0.2 ~ 21960"))
blog(format("Static series-reluctance (OFF): mu_perp = 1/(0.8/mu_x + 0.2/1)"))
blog(format("  If mu_x=1 (BH placeholder): mu_perp = 1/(0.8+0.2) = 1  *** WRONG ***"))
blog(format("  If mu_x=mu_init=36600:       mu_perp = 1/(0.8/36600+0.2) ~ 4.9"))
blog("")

Pt_on,  Ps_on,  F_on  = analyze_our("OUR-ON",  1, tmp_on)
Pt_off, Ps_off, F_off = analyze_our("OUR-OFF", 0, tmp_off)
Pt_leg, Ps_leg, F_leg = analyze_case("LEGACY",    legacy)

blog("=== Summary ===")
blog(format("P_total:  ON=%.4f W  OFF=%.4f W  Legacy=%.4f W", Pt_on, Pt_off, Pt_leg))
blog(format("P_sides:  ON=%.4f W  OFF=%.4f W  Legacy=%.4f W", Ps_on, Ps_off, Ps_leg))
blog(format("|flux|:   ON=%.3e   OFF=%.3e   Legacy=%.3e Wb-t", F_on, F_off, F_leg))
blog(format("OFF vs Legacy: dP_tot=%.4f W (%.2f %%)", Pt_off-Pt_leg, 100*(Pt_off-Pt_leg)/Pt_leg))
blog(format("ON  vs Legacy: dP_tot=%.4f W (%.2f %%)", Pt_on-Pt_leg,  100*(Pt_on-Pt_leg)/Pt_leg))
blog(format("ON  vs OFF:    dP_tot=%.4f W (%.2f %%)", Pt_on-Pt_off,  100*(Pt_on-Pt_off)/Pt_off))

blog("")
blog("DIAGNOSIS:")
if abs((Pt_off-Pt_leg)/Pt_leg) > 0.05 then
  blog("  *** OFF differs from legacy by >5% -> static series-reluctance uses mu_x=1 placeholder")
  blog("      (BH initial permeability not applied in static Mu[k] path for LamType)")
else
  blog("  OFF matches legacy within 5% -> static path is correct")
end

fh = openfile(outtxt, "w")
write(fh, buf)
closefile(fh)
