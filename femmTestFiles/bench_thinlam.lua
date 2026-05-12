-- bench_thinlam.lua  (Lua 4.0)
-- Compare mo_blockintegral(31) vs FEMM-native AC eddy solve in a controlled
-- uniform-B benchmark. Geometry built programmatically.

results = ""
function save() writeto("D:\\FEMM Source\\femmTestFiles\\bench_thinlam_results.txt"); write(results); writeto() end
function out(s) results = results .. s .. "\n"; save() end

-- parameters
B0       = 1.0           -- T peak
FREQ     = 100.0         -- Hz (skin depth ~50mm at sigma=1MS/m)
SIG_MSm  = 1.0
SIG_Sm   = SIG_MSm * 1e6
T_LAM_mm = 2.0           -- slab thickness = lam thickness
W_mm     = 20.0          -- slab width
DOMAIN   = 100.0
DEPTH_mm = 1.0

t_m  = T_LAM_mm * 1e-3
W_m  = W_mm * 1e-3
D_m  = DEPTH_mm * 1e-3
Vol  = W_m * t_m * D_m
PI   = 3.14159265358979
omega = 2*PI*FREQ
P_an  = SIG_Sm * omega * omega * B0 * B0 * t_m * t_m / 24.0 * Vol
delta = sqrt(2/(SIG_Sm*4*PI*1e-7*omega))

out("================================================================")
out("Thin-lam benchmark - mo_blockintegral(31) vs native AC solve")
out("================================================================")
out(format("B0=%.3f T  f=%g Hz  sigma=%.3f MS/m", B0, FREQ, SIG_MSm))
out(format("slab: W=%g mm  t=%g mm  depth=%g mm  V=%.3e m^3",
            W_mm, T_LAM_mm, DEPTH_mm, Vol))
out(format("skin depth (mu_r=1): %.4f m  -> t/delta=%.3e (must be << 1)",
            delta, t_m/delta))
out(format("CLOSED FORM:  P_an = sigma*w^2*B^2*t^2/24 * V = %.6e W", P_an))
out("")

newdocument(0)
mi_probdef(FREQ, "millimeters", "planar", 1e-8, DEPTH_mm, 30)

mi_addmaterial("Air",       1, 1, 0, 0, 0,       0,        0, 1, 0)
mi_addmaterial("Slab_bulk", 1, 1, 0, 0, SIG_MSm, 0,        0, 1, 0)
mi_addmaterial("Slab_lam",  1, 1, 0, 0, SIG_MSm, T_LAM_mm, 0, 1, 2)

-- Prescribed A: A_z = A0 + A1*x_mm + A2*y_mm (units: Wb/m, x,y in mm).
-- For uniform B_x = B0 [T]: A_z = B0 * y_meters = B0 * y_mm * LengthConv.
-- LengthConv for "millimeters" = 1e-3, so A2 = B0 * 1e-3.
mi_addboundprop("uniformBx", 0, 0, B0*1e-3, 0, 0, 0, 0, 0, 0)

h = DOMAIN/2
mi_addnode(-h, -h); mi_addnode(h, -h); mi_addnode(h, h); mi_addnode(-h, h)
mi_addsegment(-h, -h,  h, -h)
mi_addsegment( h, -h,  h,  h)
mi_addsegment( h,  h, -h,  h)
mi_addsegment(-h,  h, -h, -h)

mi_clearselected()
mi_selectsegment(0, -h)
mi_selectsegment(h, 0)
mi_selectsegment(0, h)
mi_selectsegment(-h, 0)
mi_setsegmentprop("uniformBx", 0, 1, 0, 0)
mi_clearselected()

sw = W_mm/2
st = T_LAM_mm/2
mi_addnode(-sw, -st); mi_addnode(sw, -st); mi_addnode(sw, st); mi_addnode(-sw, st)
mi_addsegment(-sw, -st,  sw, -st)
mi_addsegment( sw, -st,  sw,  st)
mi_addsegment( sw,  st, -sw,  st)
mi_addsegment(-sw,  st, -sw, -st)

mi_addblocklabel(0, 30)
mi_clearselected()
mi_selectlabel(0, 30)
mi_setblockprop("Air", 1, 0, "<None>", 0, 0, 0)
mi_clearselected()

mi_addblocklabel(0, 0)
mi_clearselected()
mi_selectlabel(0, 0)
-- block label args: (name, automesh, meshsize, circuit, magdir, group, turns)
mi_setblockprop("Slab_bulk", 0, 0.3, "<None>", 0, 1, 0)
mi_clearselected()

mi_saveas("D:\\FEMM Source\\femmTestFiles\\bench_thinlam.fem")

out("---- PASS A: bulk conductor (FEMM-native AC eddy solve) ----")
out("step: analyze A")
mi_analyze(1)
out("step: loadsolution A")
mi_loadsolution()

mo_clearblock()
mo_groupselectblock(1)
P_ref_resistive = abs(mo_blockintegral(4))
P_ref_hyst      = abs(mo_blockintegral(3))
Vol_femm        = abs(mo_blockintegral(10))
mo_clearblock()
mo_close()

out(format("  P(4) resistive  = %.6e W", P_ref_resistive))
out(format("  P(3) hyst+lam-mu= %.6e W", P_ref_hyst))
out(format("  Volume(10)      = %.6e m^3 (closed form V = %.6e)", Vol_femm, Vol))

out("---- PASS B: LamType=2 + mo_blockintegral(31) ----")
mi_clearselected()
mi_selectlabel(0, 0)
mi_setblockprop("Slab_lam", 0, 0.3, "<None>", 0, 1, 0)
mi_clearselected()
mi_saveas("D:\\FEMM Source\\femmTestFiles\\bench_thinlam.fem")
out("step: analyze B")
mi_analyze(1)
out("step: loadsolution B")
mi_loadsolution()

mo_clearblock()
mo_groupselectblock(1)
P_test_31 = abs(mo_blockintegral(31))
P_test_3  = abs(mo_blockintegral(3))
P_test_4  = abs(mo_blockintegral(4))
mo_clearblock()
mo_close()

out(format("  P(31) thin-lam  = %.6e W", P_test_31))
out(format("  P(3)  hyst+lam  = %.6e W", P_test_3))
out(format("  P(4)  resistive = %.6e W", P_test_4))

out("=== COMPARISON ===")
out(format("  P_an  = %.6e W", P_an))
out(format("  P_ref = %.6e W", P_ref_resistive))
out(format("  P_31  = %.6e W", P_test_31))
out(format("  P_ref / P_an  = %.4f", P_ref_resistive / P_an))
out(format("  P_31  / P_an  = %.4f", P_test_31 / P_an))
if P_ref_resistive > 0 then
    out(format("  P_31  / P_ref = %.4f", P_test_31 / P_ref_resistive))
end

save()
exit()
