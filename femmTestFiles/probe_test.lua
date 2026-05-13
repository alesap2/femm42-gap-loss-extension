src = "D:\FEMM Source\femmTestFiles\pourleroi_cc_magnetostatic_rev3.fem"
tmp = "D:\FEMM Source\femmTestFiles\probe_test_tmp.fem"
open(src)
mi_saveas(tmp)
mi_probdef(1000, "millimeters", "planar", 1e-008, 35, 30)
mi_analyse(0)
mi_loadsolution()
A,B1,B2 = mo_getpointvalues(7, 30)
fp = openfile("D:\FEMM Source\femmTestFiles\probe_test_out.txt","w")
write(fp, format("A=%.6e  Bx=%.6e  By=%.6e\n", abs(A), abs(B1), abs(B2)))
closefile(fp)
mo_close()
mi_close()
quit()
