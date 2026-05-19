-- validate_reader.lua  (Lua 4 / FEMM 4.2 compatible)
-- Loads _gap_battery4_tmp solution and reports blockintegral values
-- plus mu values at sample points.

OUT_PATH = "D:/FEMM Source/femmTestFiles/validate_reader_out.txt"

function wlog(s)
  local fp = openfile(OUT_PATH, "a")
  if fp ~= nil then write(fp, s, "\n") ; closefile(fp) end
end

-- Clear output file
local fp0 = openfile(OUT_PATH, "w")
if fp0 ~= nil then write(fp0, "") ; closefile(fp0) end

wlog("Script started")

open("D:/FEMM Source/femmTestFiles/_gap_battery4_tmp.fem")
wlog("fem opened")

mi_loadsolution()
wlog("solution loaded")

-- Sample mu at center of inner lower gap block (4.4, 17.15 mm)
local mu1x, mu1y = mo_getmu(4.4, 17.15)
wlog(format("mo_getmu(4.4,17.15): mux=%.4e  muy=%.4e", mu1x, mu1y))

-- Sample B and H at same point
local vals = mo_getpointvalues(4.4, 17.15)
-- returns: A Bx By sig e Hx Hy Je Js Bx_pu By_pu Mu_x Mu_y
wlog(format("  Bx=%.6f By=%.6f  Hx=%.4e Hy=%.4e", vals[2], vals[3], vals[6], vals[7]))

-- gap_extra = 4 mm  ->  upper label y = 53.45 + 4 = 57.45
mo_selectblock(4.4, 17.15)
mo_selectblock(4.4, 57.45)
wlog("blocks selected")

local P3  = mo_blockintegral(3)
local P4  = mo_blockintegral(4)
local P6  = mo_blockintegral(6)
local P31 = mo_blockintegral(31)
local A5  = mo_blockintegral(5)

wlog(format("blockintegral(3)  = %.6f W", P3))
wlog(format("blockintegral(4)  = %.6f W", P4))
wlog(format("blockintegral(6)  = %.6f W  [= 3+4]", P6))
wlog(format("blockintegral(31) = %.6f W  [thin-lam Bx]", P31))
wlog(format("blockintegral(5)  = %.6e m2", A5))

mo_close()
wlog("done")
