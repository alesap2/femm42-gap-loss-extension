"""Quick test: LT1 (left leg, x=7) vs LT2 (right leg, x=37) in rev5."""
import sys, subprocess, shutil
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from femm_ans_reader import FemmAns

ROOT    = Path(r"D:\FEMM Source")
FEMM    = ROOT / "TestBin" / "femm.exe"
SRC     = ROOT / "femmTestFiles" / "pourleroi_cc_magnetostatic_rev5.fem"
TMP_FEM = ROOT / "femmTestFiles" / "_rev5_test.fem"
TMP_ANS = TMP_FEM.with_suffix(".ans")
LUA     = ROOT / "femmTestFiles" / "_rev5_test.lua"

shutil.copy(SRC, TMP_FEM)

tmp_fem_lua = str(TMP_FEM).replace("\\", "\\\\")
lua_text = (
    'open("' + tmp_fem_lua + '")\n'
    'mi_probdef(100000, "millimeters", "planar", 1e-8, 35, 30)\n'
    'mi_modifycircprop("coil", 1, 300)\n'
    'mi_analyze(1)\n'
    'quit()\n'
)
LUA.write_text(lua_text, encoding="ascii")

print("Solving rev5 at 100 kHz, I=300 A...")
subprocess.run([str(FEMM), f"-lua-script={LUA}"], cwd=str(ROOT), check=True)
print("Done. Reading .ans...")

ans = FemmAns(TMP_ANS)

# Report frequency and depth
print(f"  freq={ans._freq} Hz  depth={ans._depth_m*1000:.1f} mm")

for label, (seed_x, seed_y) in [("LT1 (x=7)",  (7,  35.3)),
                                  ("LT2 (x=37)", (37, 35.3))]:
    lbl_set = ans.select_by_points([(seed_x, seed_y)])
    lbl_idx = next(iter(lbl_set))
    mat = ans._materials[ans._block_labels[lbl_idx][2]]
    Px, Py = ans.block_integral_3_split(label_indices=lbl_set)
    fbx    = ans.f_bx(label_indices=lbl_set)
    print(f"\n{label}  material={mat.name!r}  LamType={mat.lam_type}  PerpLenz={mat.perp_lenz}")
    print(f"  Px={Px*1000:.4f} mW   Py={Py*1000:.2f} mW   Ptot={( Px+Py)*1000:.2f} mW")
    print(f"  f_bx={fbx*100:.3f}%")
