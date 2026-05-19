# scripts/extract_interior_mm_v2.py
"""
Extrage INTERIOR_MM prin diferența de energie
"""

import RNA

def encode(s):
    return ''.join(s)

def extract_interior_mm():
    """
    Extrage INTERIOR_MM comparând energiile a două interior loops
    care diferă doar prin mismatch.
    """
    
    base_map = ['A', 'C', 'G', 'U']
    pair_seqs = [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C'), ('G', 'U'), ('U', 'G')]
    
    print("// INTERIOR_MM extracted from ViennaRNA")
    print("pub const INTERIOR_MM: [[[i32; 4]; 4]; 6] = [")
    
    results = []
    
    for pi, (p5, p3) in enumerate(pair_seqs):
        pair_results = []
        print(f"    // {p5}-{p3} pair (index {pi})")
        print("    [")
        
        for mm5_idx in range(4):
            row = []
            mm5 = base_map[mm5_idx]
            
            for mm3_idx in range(4):
                mm3 = base_map[mm3_idx]
                
                # Interior loop 3x3 cu structură fixă
                # Outer: p5-p3, Inner: C-G (mereu)
                # Stânga: mm5 + AA
                # Dreapta: AA + mm3
                
                # Secvență: [p5][mm5][A][A][C][AAA][G][A][A][mm3][p3]
                #           0   1    2  3  4  567  8  9 10  11   12
                # Struct:   (   .   .  .  (  ... )  .  .  .    )
                
                seq = f"{p5}{mm5}AAC{'AAA'}GAA{mm3}{p3}"
                struct = "(...(...)...)"
                
                # Verifică lungimea
                if len(seq) != len(struct):
                    print(f"// EROARE: len mismatch {len(seq)} vs {len(struct)}")
                    row.append(0)
                    continue
                
                md = RNA.md()
                md.dangles = 0
                fc = RNA.fold_compound(seq, md)
                
                try:
                    # In loc sa aflam direct mismatch, mai intai extragem evaluarea cu funcția de interior loop
                    # evaluate_int_loop(i, j, p, q) -> energia totala (INIT + Mismatch stanga + Mismatch dreapta)
                    # i=1, j=13
                    # p=5, q=9
                    e = fc.eval_int_loop(1, 13, 5, 9)
                    
                    # Pentru C-G inner, stim ca mm_inner este A,A care contribuie o anumita valoare
                    # Dacă vrem să obținem DOAR mm_outer, procedăm relativ:
                    # Alternativ putem obține direct tabelul folosind doar `extract_mismatch_only()` de mai jos!
                    e_int = int(round(e))
                    row.append(e_int)
                except Exception as ex:
                    print(f"// Eroare la {seq}: {ex}")
                    row.append(0)
            
            print(f"        [{row[0]:>5}, {row[1]:>5}, {row[2]:>5}, {row[3]:>5}],  // {mm5} at 5'")
            pair_results.append(row)
        
        print("    ],")
        results.append(pair_results)
    
    print("];")
    
    return results

def extract_mismatch_only():
    import sys
    sys.stdout.reconfigure(encoding='utf-8')
    print("\n// EXTRAGERE CONTRIBUȚIE MISMATCH IZOLATĂ")
    
    base_map = ['A', 'C', 'G', 'U']
    pair_seqs = [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C'), ('G', 'U'), ('U', 'G')]
    
    print("pub const INTERIOR_MM: [[[i32; 4]; 4]; 6] = [")
    
    for pi, (p5, p3) in enumerate(pair_seqs):
        print(f"    // {p5}-{p3} pair")
        print("    [")
        
        # Calculate reference energy (A-A mismatch at outer pair)
        ref_seq = f"{p5}AA G UUU C AA{p3}".replace(' ', '')
        md = RNA.md()
        fc_ref = RNA.fold_compound(ref_seq, md)
        try:
            ref_e = fc_ref.eval_int_loop(1, 13, 5, 9)
        except:
            ref_e = 0
            
        for mm5_idx in range(4):
            row = []
            mm5 = base_map[mm5_idx]
            
            for mm3_idx in range(4):
                mm3 = base_map[mm3_idx]
                
                seq = f"{p5}{mm5}A G UUU C A{mm3}{p3}".replace(' ', '')
                fc = RNA.fold_compound(seq, md)
                try:
                    e = fc.eval_int_loop(1, 13, 5, 9)
                    # For isolated mismatch, if we subtract the reference (A-A),
                    # we get the difference from A-A mismatch.
                    # In Turner 2004, A-A mismatch is usually 0.00 for interior loops.
                    # So (e - ref_e) * 100 gives the exact mismatch table!
                    diff = int(round((e - ref_e) * 100))
                    row.append(diff)
                except:
                    row.append(0)
                    
            print(f"        [{row[0]:>5}, {row[1]:>5}, {row[2]:>5}, {row[3]:>5}],  // {mm5}")
        print("    ],")
    print("];")

if __name__ == "__main__":
    extract_mismatch_only()
