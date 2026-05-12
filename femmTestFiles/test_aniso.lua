-- test_aniso.lua
-- Tests mi_setmataniso() and mo_blockintegral(31) for anisotropic conductivity.
-- Geometry: pourleroi_cc_magnetostatic_rev2.fem
--   "Amorphous"     LamType=0  (in-plane lams, main core yoke/columns)
--   "Amorphous gap" LamType=2  (lams || Y, gap-fringe region: only Bx transverse)
--
-- Lua 4.0 syntax (no string.format -> format, no math.sqrt -> sqrt, etc.)
-- Run from FEMM Lua console:
--   dofile("D:\\FEMM Source\\femmTestFiles\\test_aniso.lua")

fempath = "D:\\FEMM Source\\femmTestFiles\\"
outfile = fempath .. "test_aniso_results.txt"
femfile = fempath .. "pourleroi_cc_magnetostatic_rev2.fem"

results = ""
function out(s)
    print(s)
    results = results .. s .. "\n"
end
function save_results()
    writeto(outfile)
    write(results)
    writeto()
end

out("=== test_aniso.lua  " .. date() .. " ===")
out("")

-- ---------------------------------------------------------------
-- Material parameters (same bulk properties for both materials)
-- ---------------------------------------------------------------
-- sigma_m  = bulk conductivity of the amorphous metal [MS/m]
-- F        = lamination fill factor
-- t_lam    = lamination thickness [mm]
-- W_core   = cross-section width of lamination stack [mm]
--            = column width of this geometry
--
-- Homogenisation formulae:
--   sigma_t [MS/m] = F * sigma_m
--   sigma_n [S/m]  = (t_lam/W_core)^2 * sigma_m * 1e6 / F
--
sigma_m  = 0.7692308   -- MS/m (bulk amorphous)
F        = 0.80
t_lam    = 0.023       -- mm
W_core   = 14.0        -- mm  (column width of this inductor)

sigma_t = F * sigma_m
sigma_n = (t_lam / W_core)^2 * sigma_m * 1e6 / F

out(format("Material parameters:"))
out(format("  sigma_m  = %.6f MS/m", sigma_m))
out(format("  LamFill  = %.2f", F))
out(format("  t_lam    = %.3f mm", t_lam))
out(format("  W_core   = %.1f mm", W_core))
out(format("  --> sigma_t = %.5f MS/m", sigma_t))
out(format("  --> sigma_n = %.3f S/m", sigma_n))
out(format("  (for comparison: paper nanocrys. ~0.46 S/m)"))
out("")

-- ---------------------------------------------------------------
-- 1. Open pre-processor
-- ---------------------------------------------------------------
open(femfile)
mi_saveas(femfile)
out("1. Opened: " .. femfile)

-- ---------------------------------------------------------------
-- 2. Apply anisotropic conductivity
-- ---------------------------------------------------------------
out("")
out("2. Setting anisotropic conductivity...")

-- "Amorphous" (LamType=0): in-plane lams — both Bx and By drive eddies
mi_setmataniso("Amorphous", sigma_t, sigma_n)
out(format("   Amorphous     (LamType=0): sigma_t=%.5f MS/m  sigma_n=%.3f S/m", sigma_t, sigma_n))

-- "Amorphous gap" (LamType=2): lams || Y — only Bx (fringe field) is transverse
mi_setmataniso("Amorphous gap", sigma_t, sigma_n)
out(format("   Amorphous gap (LamType=2): sigma_t=%.5f MS/m  sigma_n=%.3f S/m", sigma_t, sigma_n))

-- ---------------------------------------------------------------
-- 3. Mesh & Analyze
-- ---------------------------------------------------------------
out("")
out("3. Meshing...")
mi_createmesh()
out("   Mesh created.")

out("4. Analyzing (f=1000 Hz, AC, planar)...")
mi_analyze(1)
mi_loadsolution()
out("   Solution loaded.")

-- ---------------------------------------------------------------
-- 4. Block integral(31): anisotropic eddy loss
-- ---------------------------------------------------------------
out("")
out("5. Block integral(31) — anisotropic eddy loss [W]:")

-- Left gap region  (block label at ~7, 35.3)
mo_clearblock()
mo_selectblock(7, 35.3)
P_gap_left = mo_blockintegral(31)
out(format("   Amorphous gap LEFT  (7, 35.3):  P31 = %.4f W", P_gap_left or 0))

-- Right gap region (block label at ~37, 35.3)
mo_clearblock()
mo_selectblock(37, 35.3)
P_gap_right = mo_blockintegral(31)
out(format("   Amorphous gap RIGHT (37, 35.3): P31 = %.4f W", P_gap_right or 0))

-- Total gap loss
P_gap_total = (P_gap_left or 0) + (P_gap_right or 0)
out(format("   Total gap regions: P31 = %.4f W", P_gap_total))

-- Main core (Amorphous, LamType=0) — sample in top yoke
mo_clearblock()
mo_selectblock(22, 63)
P_core = mo_blockintegral(31)
out(format("   Amorphous core (22, 63):        P31 = %.4f W", P_core or 0))

-- ---------------------------------------------------------------
-- 5. Standard eddy loss (integral 6) for comparison
-- ---------------------------------------------------------------
out("")
out("6. Standard eddy loss integral(6) — classical [W]:")

mo_clearblock()
mo_selectblock(7, 35.3)
mo_selectblock(37, 35.3)
P_std_gap = mo_blockintegral(6)
out(format("   Gap regions (std): P6 = %.4f W", P_std_gap or 0))

mo_clearblock()
mo_selectblock(22, 63)
P_std_core = mo_blockintegral(6)
out(format("   Core (std):        P6 = %.4f W", P_std_core or 0))

-- ---------------------------------------------------------------
-- 6. Summary
-- ---------------------------------------------------------------
out("")
out("=== SUMMARY ===")
out(format("  sigma_t = %.5f MS/m,  sigma_n = %.3f S/m,  W_core = %.1f mm", sigma_t, sigma_n, W_core))
out(format("  Aniso eddy loss in gap regions (P31): %.4f W", P_gap_total))
out(format("  Classic eddy loss in gap regions (P6): %.4f W", P_std_gap or 0))
if (P_std_gap or 0) > 0 then
    out(format("  Ratio P31/P6 in gap: %.4f", P_gap_total / P_std_gap))
end
out("")
out("=== END ===")

save_results()
out("Results saved to: " .. outfile)

-- Helper: accumulate lines, flush to file at end  (Lua 4.0)
results = ""
function out(s)
    print(s)
    results = results .. s .. "\n"
end
function save_results()
    writeto(outfile)
    write(results)
    writeto()
end

out("=== test_aniso.lua  " .. date() .. " ===")
out("")

-- ---------------------------------------------------------------
-- 1. Open pre-processor
-- ---------------------------------------------------------------
open(femfile)
mi_saveas(femfile)
out("1. Opened: " .. femfile)

-- ---------------------------------------------------------------
-- 2. Apply anisotropic conductivity (auto-compute)
-- ---------------------------------------------------------------
out("")
out("2. Setting anisotropic conductivity on 'Amorphous'...")
out("   d_lam=0.023 mm, LamFill=0.80, Sigma_bulk=0.769 MS/m")

mi_setmataniso("Amorphous", 0, 0)
out("   mi_setmataniso('Amorphous', 0, 0)  [auto-compute] --> OK")

-- Also test explicit values
sigma_t = 0.769 * 0.80   -- ~0.615 MS/m  in-plane
sigma_n = 500             -- S/m          normal (insulation-limited)
mi_setmataniso("Amorphous", sigma_t, sigma_n)
out(format("   mi_setmataniso('Amorphous', %.4f, %.1f)  [explicit] --> OK",
                  sigma_t, sigma_n))

-- ---------------------------------------------------------------
-- 3. Mesh
-- ---------------------------------------------------------------
out("")
out("3. Meshing...")
mi_createmesh()
out("   Mesh created.")

-- ---------------------------------------------------------------
-- 4. Analyze
-- ---------------------------------------------------------------
out("4. Analyzing (f=1000 Hz, AC)...")
mi_analyze(1)
mi_loadsolution()
out("   Solution loaded.")

-- ---------------------------------------------------------------
-- 5. Standard block integrals (baseline)
-- ---------------------------------------------------------------
out("")
out("5. Standard block integrals on 'Amorphous' core:")

mo_clearblock()
mo_selectblock(0, 30)   -- point inside the Amorphous core ring (adjust if needed)

-- mo_getpointvalues returns 14 values:
--   A, B1, B2, sigma, E, H1, H2, Je, Js, mu1, mu2, Pe, Ph, ff
A, B1, B2, sig, E, H1, H2, Je, Js, mu1, mu2, Pe, Ph, ff = mo_getpointvalues(0, 30)
Bmod = sqrt((B1 or 0)^2 + (B2 or 0)^2)
Hmod = sqrt((H1 or 0)^2 + (H2 or 0)^2)
out(format("   Point (0,30mm): |B|=%.4f T   |H|=%.2f A/m   Pe=%.4e W/m3",
                  Bmod, Hmod, Pe or 0))

A_block = mo_blockintegral(5)   -- Area [m^2]
W_eddy  = mo_blockintegral(4)   -- Ohmic (I^2R / eddy) loss [W]
W_store = mo_blockintegral(2)   -- Stored magnetic energy [J]  (NOT total loss)
out(format("   Area            = %.6e m2", A_block or 0))
out(format("   Eddy loss (std) = %.6e W   [mo_blockintegral(4)]", W_eddy  or 0))
out(format("   Stored energy   = %.6e J   [mo_blockintegral(2)]", W_store or 0))

-- ---------------------------------------------------------------
-- 6. NEW: Gap Loss Wang 2015 -- blockintegral type 31
-- ---------------------------------------------------------------
out("")
out("6. Gap Loss (anisotropic, Wang 2015):")
W_gap = mo_blockintegral(31)
out(format("   Gap Loss        = %.6e W   [mo_blockintegral(31)]", W_gap or 0))

if W_gap ~= nil and W_gap > 0 then
    out("   --> bAnisoConductivity active: PASS")
    if W_eddy ~= nil and W_eddy > 0 then
        out(format("   Ratio gap/std_eddy = %.4f", W_gap / W_eddy))
    end
else
    out("   --> W_gap = 0  (flag not set, or no AC flux in selected region)")
end

-- ---------------------------------------------------------------
-- Done
-- ---------------------------------------------------------------
out("")
out("=== TEST COMPLETE ===")
out("Results: " .. outfile)
save_results()
print("Log saved: " .. outfile)

