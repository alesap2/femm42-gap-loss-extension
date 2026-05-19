lines = open(r'D:\FEMM Source\femmTestFiles\_gap_battery4_tmp.ans').readlines()
in_gap = False
for l in lines:
    if '"Amorphous gap"' in l:
        in_gap = True
    if in_gap:
        print(l, end='')
    if in_gap and '<EndBlock>' in l:
        break
