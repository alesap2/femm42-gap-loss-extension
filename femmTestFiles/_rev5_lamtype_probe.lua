src    = "D:\\FEMM Source\\femmTestFiles\\pourleroi_cc_magnetostatic_rev5.fem"
tmpfem = "D:\\FEMM Source\\femmTestFiles\\_rev5_lamtype_probe_tmp.fem"
outtxt = "D:\\FEMM Source\\femmTestFiles\\_rev5_lamtype_probe.txt"

buf = ""
function blog(s) buf = buf .. s .. "\n" end

function readall(path)
  fp = openfile(path, "r")
  txt = ""
  line = read(fp, "*l")
  while line ~= nil do
    txt = txt .. line .. "\n"
    line = read(fp, "*l")
  end
  closefile(fp)
  return txt
end

function patch_mat_lamtype(txt, matname, lamtype)
  i1, i2 = strfind(txt, "\"" .. matname .. "\"")
  if i1 == nil then return txt end
  e1, e2 = strfind(txt, "<EndBlock>", i2)
  blk = strsub(txt, i2, e1)
  blk = gsub(blk, "<LamType> = %d+", "<LamType> = " .. lamtype, 1)
  blk = gsub(blk, "\n%s*<LamHybridSigmaZ> = %d+", "", 1)
  return strsub(txt, 1, i2-1) .. blk .. strsub(txt, e1+1)
end

function write_case(lamtype)
  txt = readall(src)
  txt = patch_mat_lamtype(txt, "Amorphous laminated", lamtype)
  txt = patch_mat_lamtype(txt, "Amorphous gap LT2", lamtype)
  fp = openfile(tmpfem, "w")
  write(fp, txt)
  closefile(fp)
end

function run_case(tag, lamtype)
  write_case(lamtype)
  open(tmpfem)
  mi_saveas(tmpfem)
  mi_analyze(1)
  mi_loadsolution()

  ci, vi, fi = mo_getcircuitproperties("coil")
  A1, Bx1, By1 = mo_getpointvalues(4.4, 17.15)
  A2, Bx2, By2 = mo_getpointvalues(22.0, 35.3)
  blog(format("%s: |flux|=%.6e |Z|=%.6e", tag, abs(fi), abs(vi)/abs(ci)))
  blog(format("  probe side gap  (4.4,17.15): |Bx|=%.6e |By|=%.6e", abs(Bx1), abs(By1)))
  blog(format("  probe centre gap (22.0,35.3): |Bx|=%.6e |By|=%.6e", abs(Bx2), abs(By2)))
  mo_close()
  mi_close()
end

function run_existing()
  open(src)
  mi_saveas(tmpfem)
  mi_analyze(1)
  mi_loadsolution()
  ci, vi, fi = mo_getcircuitproperties("coil")
  A1, Bx1, By1 = mo_getpointvalues(4.4, 17.15)
  A2, Bx2, By2 = mo_getpointvalues(22.0, 35.3)
  blog(format("as-is mixed file: |flux|=%.6e |Z|=%.6e", abs(fi), abs(vi)/abs(ci)))
  blog(format("  probe side gap  (4.4,17.15): |Bx|=%.6e |By|=%.6e", abs(Bx1), abs(By1)))
  blog(format("  probe centre gap (22.0,35.3): |Bx|=%.6e |By|=%.6e", abs(Bx2), abs(By2)))
  mo_close()
  mi_close()
end

blog("rev5 lamtype probe, hybrid sigma_z disabled")
run_existing()
run_case("all LamType=1", 1)
run_case("all LamType=2", 2)

fp = openfile(outtxt, "w")
write(fp, buf)
closefile(fp)
quit()
