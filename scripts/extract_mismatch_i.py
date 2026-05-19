import RNA

PAIRS = [('A','U'), ('U','A'), ('C','G'), ('G','C'), ('G','U'), ('U','G')]
BASES = ['A', 'C', 'G', 'U']

# INTERIOR_INIT[6] = 200 (2.0 kcal/mol)
# Asymmetry = 0
# Un 3x3 loop are structura: ((...(...)...))
# fc.eval_int_loop() = 200 + Mismatch(o) + Mismatch(i)
# Vrem să extragem Mismatch(o).
# Dacă fixăm Mismatch(i) să fie C-G cu mm 0, îl putem scădea.

def eval_3x3(oi_idx, m5, m3):
    o = PAIRS[oi_idx]
    # Folosim C-G ca inner pair (index 2)
    # Folosim A si A ca mismatch-uri inner
    # G + o[0] + m5 + A + A + C + G + A + A + m3 + o[1] + C
    # 0   1      2    3   4   5   6   7   8   9    10     11
    # Outer pair: 1, 10
    # Inner pair: 5, 6
    seq = f"G{o[0]}{BASES[m5]}AACGAA{BASES[m3]}{o[1]}C"
    fc = RNA.fold_compound(seq)
    return fc.eval_int_loop(2, 11, 6, 7) # Indexare 1-based pt ViennaRNA

# Gasim valoarea de baza (INIT + inner mismatch) folosind o valoare pe care o aflam prin diferenta.
# Din Turner 2004 stim ca mismatch_i pentru A-U cu mm A, A este un anumit numar.
# Cel mai simplu: aflam valorile relative si le centram.
# Sau chiar mai simplu, folosim un mic loop 1xn unde doar un mismatch se aplica?
# In ViennaRNA, mismatch-urile se aplica pt 1x3? Nu, 1x(n>2) primeste doar un singur mismatch! (cel de la capatul opus lui 1 nu are?).
# De fapt, formula are Mismatch_5' si Mismatch_3'.

# Solutia cea mai curată: in ViennaRNA, mismatch parameters pentru un pair X cu baze Y,Z
# este o tabela independenta. 
import sys
# Daca nu o putem extrage perfect absolut, o extragem relativ la C-G cu A,A si o ajustam ca INTERIOR_MM[C-G][A][A] = 0.
# Asta ar functiona, dar poate muta o constanta in INTERIOR_INIT.
# Dar noi am copiat INTERIOR_INIT = 200 din params.rs (care e Turner 2004).

base_val = eval_3x3(2, 0, 0) # C-G cu mm A,A 
# Stiu ca din tabelul Turner, C-G cu A,A pt interior loop are o anumita energie.
# Daca presupunem ca HAIRPIN_MM e oarecum similar, HAIRPIN_MM[C-G][A][A] = -150.
# Haide pur si simplu sa scriem un C file pe care-l compila python prin ctypes ca sa acceseze `mismatchI`? E prea complicat.
# Vienna are un fisier rna_turner2004.par în site-packages/RNA ?

            # Construim secventa G + cp_o + m5 + A + m3_in + cp_i + G + C + cp_i + m5_in + A + m3 + cp_o + C
            # Mai simplu: evaluam loop-ul cu RNA.fold_compound si determinam valoarea exactă
            pass

# Alternativă: putem descărca mismatch_i.txt din distribuția ViennaRNA!
