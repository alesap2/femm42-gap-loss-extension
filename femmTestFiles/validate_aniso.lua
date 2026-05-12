-- validate_aniso.lua  (Lua 4.0)  rev4
-- Validates anisotropic ELECTRICAL CONDUCTIVITY (Wang 2017) for LamType=2.
-- Expected (Wang eq. 8):  P_gap ∝ sigma_t · f^~1.72..2 · B^2

fempath = "D:\\FEMM Source\\femmTestFiles\\"
outfile = fempath .. "validate_aniso_results.txt"
femfile = fempath .. "pourleroi_cc_magnetostatic_rev3.fem"

results = ""
function out(s) print(s); results = results .. s .. "\n" end
function save() writeto(outfile); write(results); writeto() end

-- The "Amorphous gap" block labels in rev3.fem (BlockType=3 in the .fem file)
gap_pts_x = { 4.4, 4.4, 38.7, 37.6 }
gap_pts_y = {17.15, 58.45, 23.3, 50.7 }
N_GAP = 4

-- Copper coil labels (LitzCu_0p2) — 16 turns
cu_x = { 26.425, 47.575, 26.425, 47.575, 26.425, 47.575, 26.425, 47.575,
         17.575, -3.575, 17.575, -3.575, 17.575, -3.575, 17.575, -3.575 }
cu_y = { 19.325, 19.325, 26.975, 26.975, 48.625, 48.625, 56.275, 56.275,
         19.325, 19.325, 26.975, 26.975, 48.625, 48.625, 56.275, 56.275 }

-- Amorphous core labels (BH curve, sigma=0.7692, LamType=2 in-plane)
core_x = { 4.4, 4.4, 21.2, 38.7, 37.6, 21.9,
           231.21411, 221.87078, 204.39474, 179.25348, 147.23426, 109.42075, 67.157 }
core_y = { 17.15, 58.45, 71.3, 23.3, 50.7, 3.8,
           76.91527, 118.08919, 157.17227, 192.55348, 222.72632, 246.35236, 262.31957 }

function select_pts(xs, ys, n)
    mo_clearblock()
    local kk = 1
    while kk <= n do
        mo_selectblock(xs[kk], ys[kk])
        kk = kk + 1
    end
end

function select_gap()  select_pts(gap_pts_x, gap_pts_y, N_GAP)  end
function select_cu()   select_pts(cu_x, cu_y, 16) end
function select_core() select_pts(core_x, core_y, 13) end

-- BlockIntegral types:
--   2  = total losses
--   4  = total volume (sanity)
--   6  = volume of selected blocks
--  17  = real-power dissipation (Re part of A·J*)
--  31  = thin-lam analytical formula  σ_t·ω²·B_perp²·t²/24·F · Vol
--  (FEMM types: 0=A·J real, 1=imag, 2=total losses, etc. See manual)

function probe_run(label, freq, sigma_t_MSm, sigma_n_Sm)
    open(femfile)
    mi_probdef(freq)                              -- set frequency
    mi_setmataniso("Amorphous gap", sigma_t_MSm, sigma_n_Sm)
    mi_saveas(fempath .. "v_tmp.fem")
    mi_createmesh()
    mi_analyze(1)
    mi_loadsolution()

    select_gap()
    P_gap_h   = abs(mo_blockintegral(3))   -- hysteresis + laminated eddy (mu_im)
    P_gap_r   = abs(mo_blockintegral(4))   -- macro resistive (0 for laminated)
    P_gap_th  = abs(mo_blockintegral(31))  -- our explicit thin-lam term
    P_thin31_keep = P_gap_th
    Vol_gap   = abs(mo_blockintegral(10))

    select_core()
    P_core_h  = abs(mo_blockintegral(3))
    P_core_r  = abs(mo_blockintegral(4))

    select_cu()
    P_cu      = abs(mo_blockintegral(4))

    I, V, Flx = mo_getcircuitproperties("coil")
    Z = V/I
    Rh = re(Z); Xh = im(Z)
    Lh = Xh / (2*3.14159265358979*freq)
    P_circ = 0.5 * Rh * abs(I) * abs(I)

    out(format("  %-32s f=%6.0fHz sig_t=%.3eMS/m sig_n=%.3eS/m", label, freq, sigma_t_MSm, sigma_n_Sm))
    out(format("    P_gap_h(3)=%.4e  P_gap_r(4)=%.4e  P_thin31=%.4e", P_gap_h, P_gap_r, P_gap_th))
    out(format("    P_core_h(3)=%.4e P_core_r(4)=%.4e P_cu(4)=%.4e P_RI2/2=%.4e",
                P_core_h, P_core_r, P_cu, P_circ))
    out(format("    R=%.4e Ohm   L=%.4e uH   |I|=%.3e   Vol_gap=%.4e",
                Rh, Lh*1e6, abs(I), Vol_gap))

    mo_close(); mi_close()
    return P_thin31_keep, P_gap_h, P_core_h, P_cu
end

out("=== validate_aniso.lua " .. date() .. " ===")
out("Validating ANISOTROPIC CONDUCTIVITY (Wang 2017 model)")
out("Base material: sigma=0.7692 MS/m, t_lam=0.023mm, F=0.8, LamType=2")
out("")

-- ============================================================
-- TEST 1: sigma_t linear scaling at f=1kHz
-- ============================================================
out("--- TEST 1: P_gap vs sigma_t (linear scaling expected) ---")
S1_05 = probe_run("T1: sigma_t x0.5",  1000, 0.3077, 1.0)
S1_10 = probe_run("T1: sigma_t x1.0 (base)", 1000, 0.6154, 1.0)
S1_20 = probe_run("T1: sigma_t x2.0",  1000, 1.2308, 1.0)
S1_40 = probe_run("T1: sigma_t x4.0",  1000, 2.4616, 1.0)
out("")
out(format("  Ratio (x2.0)/(x1.0) = %.3f   expected ~2.0", S1_20/S1_10))
out(format("  Ratio (x4.0)/(x1.0) = %.3f   expected ~4.0", S1_40/S1_10))
out(format("  Ratio (x0.5)/(x1.0) = %.3f   expected ~0.5", S1_05/S1_10))
out("")

-- ============================================================
-- TEST 2: frequency scaling at fixed sigma_t
-- ============================================================
out("--- TEST 2: P_gap vs frequency (alpha ~ 2.0 thin-lam, ->1.72 high f) ---")
F_list = {100, 300, 1000, 3000, 10000, 30000, 100000}
P_at_f = {}
k = 1
while k <= 7 do
    f = F_list[k]
    P_at_f[k] = probe_run(format("T2: f=%dHz", f), f, 0.6154, 1.0)
    k = k + 1
end
out("")
out("  Fitted exponents (log-log slope P_gap vs f):")
k = 1
while k < 7 do
    if P_at_f[k] == nil or P_at_f[k+1] == nil or P_at_f[k] <= 0 or P_at_f[k+1] <= 0 then
        out(format("    %6dHz -> %6dHz:  SKIP (nil or zero)", F_list[k], F_list[k+1]))
    else
        alpha = log(P_at_f[k+1]/P_at_f[k]) / log(F_list[k+1]/F_list[k])
        out(format("    %6dHz -> %6dHz:  alpha = %.3f", F_list[k], F_list[k+1], alpha))
    end
    k = k + 1
end
out("")

-- ============================================================
-- TEST 3: sigma_n influence at f=1kHz
-- ============================================================
out("--- TEST 3: P_gap vs sigma_n (small effect for thin lams expected) ---")
N1 = probe_run("T3: sigma_n=0.01 S/m", 1000, 0.6154, 0.01)
N2 = probe_run("T3: sigma_n=1 S/m   ", 1000, 0.6154, 1.0)
N3 = probe_run("T3: sigma_n=100 S/m ", 1000, 0.6154, 100.0)
N4 = probe_run("T3: sigma_n=10000   ", 1000, 0.6154, 10000.0)
out("")
out(format("  Ratio sn=1/sn=0.01 = %.3f", N2/N1))
out(format("  Ratio sn=100/sn=1  = %.3f", N3/N2))
out(format("  Ratio sn=1e4/sn=1  = %.3f", N4/N2))
out("")

-- ============================================================
-- TEST 4: bAnisoConductivity OFF vs ON
-- ============================================================
out("--- TEST 4: anisotropy flag OFF (sigma_t=0,sigma_n>0) vs ON ---")
-- OFF: sigma_t=0, sigma_n>0 → bAnisoConductivity=FALSE, falls back to bulk Cduct
open(femfile)
mi_probdef(1000)
mi_setmataniso("Amorphous gap", 0., 1.)
mi_saveas(fempath .. "v_tmp.fem")
mi_createmesh(); mi_analyze(1); mi_loadsolution()
select_gap()
A_off_gap = abs(mo_blockintegral(3))
A_off_thin = abs(mo_blockintegral(31))
out(format("  aniso OFF: P_gap_h=%.4e  P_thin31=%.4e", A_off_gap, A_off_thin))
mo_close(); mi_close()

A_on = probe_run("T4: aniso ON (sigma_t=0.6154)", 1000, 0.6154, 1.0)
out(format("  P_gap OFF / ON = %.3f  (P_thin31 ratio; expected = 1/LamFill = %.3f)",
            A_off_thin/A_on, 1.0/0.8))
out("")

-- ============================================================
-- TEST 5: Lamination thickness scaling
-- Wang eq. 8: P ∝ t_lam^2 (in thin-lam regime)
-- ============================================================
out("--- TEST 5: P_thin31 vs Lam_d (expect quadratic scaling t^2) ---")
function run_lam_d(t_lam_mm)
    open(femfile)
    mi_probdef(1000)
    mi_modifymaterial("Amorphous gap", 6, t_lam_mm)   -- propnum 6 = Lam_d [mm]
    mi_setmataniso("Amorphous gap", 0.6154, 1.0)
    mi_saveas(fempath .. "v_tmp.fem")
    mi_createmesh(); mi_analyze(1); mi_loadsolution()
    select_gap()
    Pt = abs(mo_blockintegral(31))
    Ph = abs(mo_blockintegral(3))
    Iv, Vv, Fv = mo_getcircuitproperties("coil")
    Lv = im(Vv/Iv) / (2*3.14159265358979*1000)
    out(format("    t_lam=%.4f mm   P_thin31=%.4e   P_gap_h(3)=%.4e   L=%.3f uH",
                t_lam_mm, Pt, Ph, Lv*1e6))
    mo_close(); mi_close()
    return Pt
end

T_list = {0.005, 0.010, 0.023, 0.050, 0.100, 0.200}
P_t = {}
ii = 1
while ii <= 6 do
    P_t[ii] = run_lam_d(T_list[ii])
    ii = ii + 1
end
out("")
out("  Fitted exponents (log-log slope P vs t_lam):")
ii = 1
while ii < 6 do
    if P_t[ii] > 0 and P_t[ii+1] > 0 then
        beta = log(P_t[ii+1]/P_t[ii]) / log(T_list[ii+1]/T_list[ii])
        out(format("    %.4f -> %.4f mm:  beta = %.3f", T_list[ii], T_list[ii+1], beta))
    end
    ii = ii + 1
end
out("")

out("=== END ===")
save()
out("Saved: " .. outfile)
