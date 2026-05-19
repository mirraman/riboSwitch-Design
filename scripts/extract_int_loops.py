import RNA
import sys

# Indexul perechilor în Rust:
# 0: A-U, 1: U-A, 2: C-G, 3: G-C, 4: G-U, 5: U-G
PAIRS = [('A','U'), ('U','A'), ('C','G'), ('G','C'), ('G','U'), ('U','G')]
BASES = ['A', 'C', 'G', 'U']

def build_seq_1x1(oi_idx, ii_idx, mm5_idx, mm3_idx):
    # i, j = 1, 6
    # p, q = 3, 4
    # G + cp_out[0] + mm5 + cp_in[0] + cp_in[1] + mm3 + cp_out[1] + C
    o = PAIRS[oi_idx]
    i = PAIRS[ii_idx]
    mm5 = BASES[mm5_idx]
    mm3 = BASES[mm3_idx]
    
    # Adaugăm un G-C la capete ca să nu fie lonely pair
    seq = f"G{o[0]}{mm5}{i[0]}AAAA{i[1]}{mm3}{o[1]}C"
    # len = 1+1+1+1 + 4 + 1+1+1+1 = 12
    # 1: G, 12: C
    # 2: o[0], 11: o[1]
    # 4: i[0], 9: i[1]
    # i=2, j=11
    # p=4, q=9
    # deci mm5 e la 3, mm3 e la 10
    return seq, 2, 11, 4, 9

def build_seq_1x2(oi_idx, ii_idx, mm5_idx, mm3_idx, mm_mid_idx):
    o = PAIRS[oi_idx]
    i = PAIRS[ii_idx]
    mm5 = BASES[mm5_idx]
    mm3 = BASES[mm3_idx]
    mm_mid = BASES[mm_mid_idx]
    
    # 1x2 are 1 pe stânga, 2 pe dreapta
    seq = f"G{o[0]}{mm5}{i[0]}AAAA{i[1]}{mm_mid}{mm3}{o[1]}C"
    # len = 1+1+1+1 + 4 + 1+1+1+1+1 = 13
    # 2: o[0], 12: o[1]
    # 4: i[0], 9: i[1]
    return seq, 2, 12, 4, 9

def build_seq_2x1(oi_idx, ii_idx, mm5_idx, mm3_idx, mm_mid_idx):
    o = PAIRS[oi_idx]
    i = PAIRS[ii_idx]
    mm5 = BASES[mm5_idx]
    mm3 = BASES[mm3_idx]
    mm_mid = BASES[mm_mid_idx]
    
    # 2x1 are 2 pe stânga, 1 pe dreapta
    seq = f"G{o[0]}{mm5}{mm_mid}{i[0]}AAAA{i[1]}{mm3}{o[1]}C"
    # len = 1+1+1+1+1 + 4 + 1+1+1+1 = 13
    # 2: o[0], 12: o[1]
    # 5: i[0], 10: i[1]
    return seq, 2, 12, 5, 10

def build_seq_2x2(oi_idx, ii_idx, mm5_out, mm3_out, mm5_in, mm3_in):

    o = PAIRS[oi_idx]
    i = PAIRS[ii_idx]
    b5o = BASES[mm5_out]
    b3o = BASES[mm3_out]
    b5i = BASES[mm5_in]
    b3i = BASES[mm3_in]
    
    seq = f"G{o[0]}{b5o}{b5i}{i[0]}AAAA{i[1]}{b3i}{b3o}{o[1]}C"
    # 2: o[0], 14: o[1]
    # 5: i[0], 10: i[1]
    return seq, 2, 14, 5, 10

INF = 10000000

print("pub const INT11: [[[[i32; 4]; 4]; 6]; 6] = [")
for oi in range(6):
    print("    [")
    for ii in range(6):
        print("        [")
        for m5 in range(4):
            vals = []
            for m3 in range(4):
                seq, a, b, p, q = build_seq_1x1(oi, ii, m5, m3)
                fc = RNA.fold_compound(seq)
                val = fc.eval_int_loop(a, b, p, q)
                if val >= 1000:
                    vals.append("INF")
                else:
                    vals.append(str(val))
            print("            [" + ", ".join(vals) + "],")
        print("        ],")
    print("    ],")
print("];")

print("pub const INT21: [[[[[i32; 4]; 4]; 4]; 6]; 6] = [")
for oi in range(6):
    print("    [")
    for ii in range(6):
        print("        [")
        for m5 in range(4):
            print("            [")
            for m3 in range(4):
                vals = []
                for mm in range(4):
                    seq, a, b, p, q = build_seq_1x2(oi, ii, m5, m3, mm)
                    fc = RNA.fold_compound(seq)
                    val = fc.eval_int_loop(a, b, p, q)
                    if val >= 1000:
                        vals.append("INF")
                    else:
                        vals.append(str(val))
                print("                [" + ", ".join(vals) + "],")
            print("            ],")
        print("        ],")
    print("    ],")
print("];")

print("pub const INT12: [[[[[i32; 4]; 4]; 4]; 6]; 6] = [")
for oi in range(6):
    print("    [")
    for ii in range(6):
        print("        [")
        for m5 in range(4):
            print("            [")
            for m3 in range(4):
                vals = []
                for mm in range(4):
                    seq, a, b, p, q = build_seq_2x1(oi, ii, m5, m3, mm)
                    fc = RNA.fold_compound(seq)
                    val = fc.eval_int_loop(a, b, p, q)
                    if val >= 1000:
                        vals.append("INF")
                    else:
                        vals.append(str(val))
                print("                [" + ", ".join(vals) + "],")
            print("            ],")
        print("        ],")
    print("    ],")
print("];")

print("pub const INT22: [[[[[[i32; 4]; 4]; 4]; 4]; 6]; 6] = [")
for oi in range(6):
    print("    [")
    for ii in range(6):
        print("        [")
        for m5o in range(4):
            print("            [")
            for m3o in range(4):
                print("                [")
                for m5i in range(4):
                    vals = []
                    for m3i in range(4):
                        seq, a, b, p, q = build_seq_2x2(oi, ii, m5o, m3o, m5i, m3i)
                        fc = RNA.fold_compound(seq)
                        val = fc.eval_int_loop(a, b, p, q)
                        if val >= 1000:
                            vals.append("INF")
                        else:
                            vals.append(str(val))
                    print("                    [" + ", ".join(vals) + "],")
                print("                ],")
            print("            ],")
        print("        ],")
    print("    ],")
print("];")
