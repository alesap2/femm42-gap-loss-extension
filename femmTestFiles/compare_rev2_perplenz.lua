-- compare_rev2_perplenz.lua  (Lua 4.0 / FEMM 4.2)
src     = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev2.fem"
tmp_on  = "D:\\FEMM Source\\femmTestFiles\\_rev2_on.fem"
tmp_off = "D:\\FEMM Source\\femmTestFiles\\_rev2_off.fem"
outtxt  = "D:\\FEMM Source\\femmTestFiles\\compare_rev2_results.txt"

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

  maxB = 0
  xmin = -10; xmax = 45; ymin = 0; ymax = 76
  nx = 50; ny = 70
  i = 0
  while i < nx do
    j = 0
    while j < ny do
      x = xmin + (xmax-xmin)*i/(nx-1)
      y = ymin + (ymax-ymin)*j/(ny-1)
      A, B1, B2 = mo_getpointvalues(x,y)
      if B1 ~= nil and B2 ~= nil then
        b1 = abs(B1); b2 = abs(B2)
        b  = sqrt(b1*b1 + b2*b2)
        if b > maxB then maxB = b end
      end
      j = j + 1
    end
    i = i + 1
  end

  ci, vi, fi = mo_getcircuitproperties("coil")
  blog(format("[%s] max|B| in window = %.4e T", tag, maxB))
  blog(format("[%s] |I_coil|         = %.4e A", tag, abs(ci)))
  blog(format("[%s] |V_coil|         = %.4e V", tag, abs(vi)))
  blog(format("[%s] |Z_coil|         = %.4e Ohm", tag, abs(vi)/abs(ci)))
  blog(format("[%s] |flux linkage|   = %.4e Wb-t", tag, abs(fi)))
  blog("")

  mo_close()
  mi_close()
  return maxB
end

blog("=== rev2 PerpLenz comparison ===")
m_on  = run_case("ON",  1, tmp_on)
m_off = run_case("OFF", 0, tmp_off)
blog(format("Ratio maxB ON/OFF = %.3f", m_on/m_off))

fh = openfile(outtxt, "w")
write(fh, buf)
closefile(fh)
