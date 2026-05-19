use crate::params::*;
const INF64: i64 = 10_000_000_000;
pub fn fold_mfe(seq: &[u8]) -> (i32, String) {
    let (mfe, _pairs, db) = fold_mfe_full(seq);
    (mfe, db)
}
pub fn fold_mfe_full(seq: &[u8]) -> (i32, Vec<i32>, String) {
    let n = seq.len();
    if n < MIN_HAIRPIN + 2 {
        return (0, vec![-1; n], ".".repeat(n));
    }
    let mut v  = vec![INF64; n * n];
    let mut wm = vec![INF64; n * n];
    macro_rules! V  { ($i:expr,$j:expr) => { v [$i*n+$j] } }
    macro_rules! WM { ($i:expr,$j:expr) => { wm[$i*n+$j] } }
    for span in (MIN_HAIRPIN + 2)..=n {
        for i in 0..=(n - span) {
            let j = i + span - 1;
            if can_pair(seq, i, j) {
                V!(i, j) = fill_v(seq, &v, &wm, i, j, n);
            }
            WM!(i, j) = fill_wm(seq, &v, &wm, i, j, n);
        }
    }
    let f5 = fill_f5(seq, &v, n);
    let mfe = f5[n] as i32;
    let mut pairs = vec![-1i32; n];
    trace_f5(seq, &v, &wm, &f5, n, &mut pairs);
    let mut db = vec![b'.'; n];
    for i in 0..n {
        if pairs[i] > i as i32 {
            db[i] = b'(';
            db[pairs[i] as usize] = b')';
        }
    }
    (mfe, pairs, String::from_utf8(db).unwrap())
}
fn fill_f5(seq: &[u8], v: &[i64], n: usize) -> Vec<i64> {
    let mut f5 = vec![0i64; n + 1];
    for i in 1..=n {
        let mut best = f5[i - 1];
        for k in 0..i {
            if i - 1 - k < MIN_HAIRPIN + 1 { continue; }
            let v_ki = v[k * n + (i - 1)];
            if v_ki >= INF64 { continue; }
            let pi = match pair_index(seq[k], seq[i - 1]) {
                Some(x) => x, None => continue,
            };
            let mut e = f5[k] + v_ki;
            if is_au_gu(pi) { e += TERMINAL_AU_PENALTY as i64; }
            if k > 0 { e += DANGLE5[pi][seq[k - 1] as usize] as i64; }
            if i < n { e += DANGLE3[pi][seq[i] as usize] as i64; }
            if e < best { best = e; }
        }
        f5[i] = best;
    }
    f5
}
fn fill_v(seq: &[u8], v: &[i64], wm: &[i64], i: usize, j: usize, n: usize) -> i64 {
    let pi = match pair_index(seq[i], seq[j]) { Some(x) => x, None => return INF64 };
    let mut best = hairpin64(seq, i, j, pi);
    let max_int = (j - i - 2).min(30);
    for p in i + 1..j {
        let nl = p - i - 1;
        if nl > max_int { break; }
        for q in (p + 1..j).rev() {
            let nr = j - q - 1;
            if nl + nr > max_int { continue; }
            if !can_pair(seq, p, q) { continue; }
            if q - p - 1 < MIN_HAIRPIN { continue; }
            let vp = v[p * n + q];
            if vp >= INF64 { continue; }
            let pp = match pair_index(seq[p], seq[q]) { Some(x) => x, None => continue };
            let e = if nl == 0 && nr == 0 {
                STACK[pi][pp] as i64 + vp
            } else if nl == 0 || nr == 0 {
                bulge64(pi, pp, nl + nr) + vp
            } else {
                interior64(seq, i, j, p, q, pi, pp, nl, nr) + vp
            };
            if e < best { best = e; }
        }
    }
    if j - i - 1 >= 2 * (MIN_HAIRPIN + 2) {
        let mut ml_base = (ML_OFFSET + ML_PER_BRANCH) as i64;
        if is_au_gu(pi) { ml_base += TERMINAL_AU_PENALTY as i64; }
        ml_base += DANGLE3[pi][seq[i + 1] as usize] as i64;
        ml_base += DANGLE5[pi][seq[j - 1] as usize] as i64;
        for k in (i + 2 + MIN_HAIRPIN)..(j - MIN_HAIRPIN - 1) {
            let wml = wm[(i + 1) * n + k];
            let wmr = wm[(k + 1) * n + j - 1];
            if wml < INF64 && wmr < INF64 {
                let e = ml_base + wml + wmr;
                if e < best { best = e; }
            }
        }
    }
    best
}
fn fill_wm(seq: &[u8], v: &[i64], wm: &[i64], i: usize, j: usize, n: usize) -> i64 {
    let mut best = INF64;
    // Evaluate single stem with optional dangles
    // Stem can be at (p, q) where p is i or i+1, and q is j or j-1
    for p in i..=i+1 {
        for q in (j.saturating_sub(1)..=j).rev() {
            if p >= q || q - p < MIN_HAIRPIN + 1 { continue; }
            if can_pair(seq, p, q) {
                let vpq = v[p * n + q];
                if vpq < INF64 {
                    let pi = pair_index(seq[p], seq[q]).unwrap();
                    let mut e = vpq + ML_PER_BRANCH as i64;
                    if is_au_gu(pi) { e += TERMINAL_AU_PENALTY as i64; }
                    
                    // Add dangle 5' if p == i + 1
                    if p == i + 1 { e += DANGLE5[pi][seq[i] as usize] as i64; }
                    // Add dangle 3' if q == j - 1
                    if q == j - 1 { e += DANGLE3[pi][seq[j] as usize] as i64; }
                    
                    // Add unpaired penalties for the dangles
                    e += ((p - i) as i64 + (j - q) as i64) * ML_PER_UNPAIRED as i64;
                    
                    if e < best { best = e; }
                }
            }
        }
    }
    if i + 1 <= j {
        let e = wm[(i + 1) * n + j];
        if e < INF64 {
            let e2 = e + ML_PER_UNPAIRED as i64;
            if e2 < best { best = e2; }
        }
    }
    if j >= 1 {
        let e = wm[i * n + j - 1];
        if e < INF64 {
            let e2 = e + ML_PER_UNPAIRED as i64;
            if e2 < best { best = e2; }
        }
    }
    for k in i + 1..j {
        let el = wm[i * n + k];
        let er = wm[(k + 1) * n + j];
        if el < INF64 && er < INF64 {
            let e = el + er;
            if e < best { best = e; }
        }
    }
    best
}
fn trace_f5(
    seq: &[u8], v: &[i64], wm: &[i64], f5: &[i64],
    n: usize, pairs: &mut Vec<i32>,
) {
    let mut i = n;
    while i > 0 {
        if f5[i] == f5[i - 1] {
            i -= 1;
            continue;
        }
        let mut matched = false;
        for k in 0..i {
            if i - 1 - k < MIN_HAIRPIN + 1 { continue; }
            let v_ki = v[k * n + (i - 1)];
            if v_ki >= INF64 { continue; }
            let pi = match pair_index(seq[k], seq[i - 1]) {
                Some(x) => x, None => continue,
            };
            let mut e = f5[k] + v_ki;
            if is_au_gu(pi) { e += TERMINAL_AU_PENALTY as i64; }
            if k > 0 { e += DANGLE5[pi][seq[k - 1] as usize] as i64; }
            if i < n { e += DANGLE3[pi][seq[i] as usize] as i64; }
            if e == f5[i] {
                pairs[k] = (i - 1) as i32;
                pairs[i - 1] = k as i32;
                trace_v(seq, v, wm, k, i - 1, n, pairs);
                i = k;
                matched = true;
                break;
            }
        }
        if !matched {
            i -= 1;
        }
    }
}
fn trace_v(
    seq: &[u8], v: &[i64], wm: &[i64],
    i: usize, j: usize, n: usize,
    pairs: &mut Vec<i32>,
) {
    let target = v[i * n + j];
    let pi = match pair_index(seq[i], seq[j]) { Some(x) => x, None => return };
    if hairpin64(seq, i, j, pi) == target { return; }
    let max_int = (j - i - 2).min(30);
    for p in i + 1..j {
        let nl = p - i - 1;
        if nl > max_int { break; }
        for q in (p + 1..j).rev() {
            let nr = j - q - 1;
            if nl + nr > max_int { continue; }
            if !can_pair(seq, p, q) { continue; }
            if q - p - 1 < MIN_HAIRPIN { continue; }
            let vp = v[p * n + q];
            if vp >= INF64 { continue; }
            let pp = match pair_index(seq[p], seq[q]) { Some(x) => x, None => continue };
            let e = if nl == 0 && nr == 0 {
                STACK[pi][pp] as i64 + vp
            } else if nl == 0 || nr == 0 {
                bulge64(pi, pp, nl + nr) + vp
            } else {
                interior64(seq, i, j, p, q, pi, pp, nl, nr) + vp
            };
            if e == target {
                pairs[p] = q as i32;
                pairs[q] = p as i32;
                trace_v(seq, v, wm, p, q, n, pairs);
                return;
            }
        }
    }
    if j - i - 1 >= 2 * (MIN_HAIRPIN + 2) {
        let mut ml_base = (ML_OFFSET + ML_PER_BRANCH) as i64;
        if is_au_gu(pi) { ml_base += TERMINAL_AU_PENALTY as i64; }
        ml_base += DANGLE3[pi][seq[i + 1] as usize] as i64;
        ml_base += DANGLE5[pi][seq[j - 1] as usize] as i64;
        for k in (i + 2 + MIN_HAIRPIN)..(j - MIN_HAIRPIN - 1) {
            let wml = wm[(i + 1) * n + k];
            let wmr = wm[(k + 1) * n + j - 1];
            if wml < INF64 && wmr < INF64 && ml_base + wml + wmr == target {
                trace_wm(seq, v, wm, i + 1, k, n, pairs);
                trace_wm(seq, v, wm, k + 1, j - 1, n, pairs);
                return;
            }
        }
    }
}
fn trace_wm(
    seq: &[u8], v: &[i64], wm: &[i64],
    i: usize, j: usize, n: usize,
    pairs: &mut Vec<i32>,
) {
    if i > j { return; }
    let target = wm[i * n + j];
    if target >= INF64 { return; }
    for p in i..=i+1 {
        for q in (j.saturating_sub(1)..=j).rev() {
            if p >= q || q - p < MIN_HAIRPIN + 1 { continue; }
            if can_pair(seq, p, q) {
                let vpq = v[p * n + q];
                if vpq < INF64 {
                    let pi = pair_index(seq[p], seq[q]).unwrap();
                    let mut e = vpq + ML_PER_BRANCH as i64;
                    if is_au_gu(pi) { e += TERMINAL_AU_PENALTY as i64; }
                    
                    if p == i + 1 { e += DANGLE5[pi][seq[i] as usize] as i64; }
                    if q == j - 1 { e += DANGLE3[pi][seq[j] as usize] as i64; }
                    e += ((p - i) as i64 + (j - q) as i64) * ML_PER_UNPAIRED as i64;
                    
                    if e == target {
                        pairs[p] = q as i32;
                        pairs[q] = p as i32;
                        trace_v(seq, v, wm, p, q, n, pairs);
                        return;
                    }
                }
            }
        }
    }
    if i + 1 <= j {
        let e = wm[(i + 1) * n + j];
        if e < INF64 && e + ML_PER_UNPAIRED as i64 == target {
            trace_wm(seq, v, wm, i + 1, j, n, pairs);
            return;
        }
    }
    if j >= 1 {
        let e = wm[i * n + j - 1];
        if e < INF64 && e + ML_PER_UNPAIRED as i64 == target {
            trace_wm(seq, v, wm, i, j - 1, n, pairs);
            return;
        }
    }
    for k in i + 1..j {
        let el = wm[i * n + k];
        let er = wm[(k + 1) * n + j];
        if el < INF64 && er < INF64 && el + er == target {
            trace_wm(seq, v, wm, i, k, n, pairs);
            trace_wm(seq, v, wm, k + 1, j, n, pairs);
            return;
        }
    }
}
#[inline]
fn can_pair(seq: &[u8], i: usize, j: usize) -> bool {
    pair_index(seq[i], seq[j]).is_some()
}
#[inline]
fn hairpin64(seq: &[u8], i: usize, j: usize, pi: usize) -> i64 {
    let e = hairpin_e(seq, i, j, pi);
    if e >= INF { INF64 } else { e as i64 }
}
fn hairpin_e(seq: &[u8], i: usize, j: usize, pi: usize) -> i32 {
    let size = j - i - 1;
    if size < MIN_HAIRPIN { return INF; }
    let mut energy = if size <= 30 {
        HAIRPIN_INIT[size]
    } else {
        HAIRPIN_INIT[30]
            + (LOOP_EXTRAPOLATION_COEFF * (size as f64 / 30.0).ln() * 100.0).round() as i32
    };
    if size == 3 {
        let w = &seq[i..=j];
        let bonus = triloop_bonus(w);
        if bonus != 0 { energy = bonus; }
        if is_au_gu(pi) { energy += TERMINAL_AU_PENALTY; }
        if (i + 1..j).all(|k| seq[k] == 1) { energy += HAIRPIN_C3; }
        return energy;
    }
    if size == 4 {
        let w = &seq[i..=j];
        let bonus = tetraloop_bonus(w);
        if bonus != 0 { return bonus; }
    }
    let b5 = seq[i + 1] as usize;
    let b3 = seq[j - 1] as usize;
    energy += HAIRPIN_MM[pi][b5][b3];
    let mm = (seq[i + 1], seq[j - 1]);
    if mm == (3, 3) || mm == (2, 0) { energy += HAIRPIN_UU_GA_BONUS; }
    if mm == (2, 2) { energy += HAIRPIN_GG_BONUS; }
    if pi == 4 && i >= 2 && seq[i - 1] == 2 && seq[i - 2] == 2 { energy += HAIRPIN_SPECIAL_GU; }
    if (i + 1..j).all(|k| seq[k] == 1) {
        energy += HAIRPIN_C_SLOPE * size as i32 + HAIRPIN_C_INTERCEPT;
    }
    energy
}
#[inline]
fn bulge64(pi: usize, pp: usize, size: usize) -> i64 {
    let mut e = if size <= 30 { BULGE_INIT[size] } else {
        BULGE_INIT[30] + (LOOP_EXTRAPOLATION_COEFF * (size as f64 / 30.0).ln() * 100.0).round() as i32
    };
    if size == 1 { e += STACK[pi][pp]; }
    if is_au_gu(pi) { e += TERMINAL_AU_PENALTY; }
    if is_au_gu(pp) { e += TERMINAL_AU_PENALTY; }
    e as i64
}
#[inline]
fn interior64(
    seq: &[u8],
    i: usize, j: usize,
    p: usize, q: usize,
    pi: usize, pp: usize,
    nl: usize, nr: usize,
) -> i64 {
    let total = nl + nr;
    if nl == 1 && nr == 1 {
        let mm5 = seq[i + 1] as usize;
        let mm3 = seq[j - 1] as usize;
        let val = INT11[pi][pp][mm5][mm3];
        if val < INF { return val as i64; }
    }
    if nl == 1 && nr == 2 {
        let mm5 = seq[i + 1] as usize;
        let mm3 = seq[j - 1] as usize;
        let mm_mid = seq[j - 2] as usize;
        let val = INT21[pi][pp][mm5][mm3][mm_mid];
        if val < INF { return val as i64; }
    }
    if nl == 2 && nr == 1 {
        let mm5 = seq[i + 1] as usize;
        let mm3 = seq[j - 1] as usize;
        let mm_mid = seq[i + 2] as usize;
        let val = INT12[pi][pp][mm5][mm3][mm_mid];
        if val < INF { return val as i64; }
    }
    if nl == 2 && nr == 2 {
        let m5o = seq[i + 1] as usize;
        let m3o = seq[j - 1] as usize;
        let m5i = seq[i + 2] as usize;
        let m3i = seq[j - 2] as usize;
        let val = INT22[pi][pp][m5o][m3o][m5i][m3i];
        if val < INF { return val as i64; }
    }
    let mut e = if total <= 30 { INTERIOR_INIT[total] } else {
        INTERIOR_INIT[30] + (LOOP_EXTRAPOLATION_COEFF * (total as f64 / 30.0).ln() * 100.0).round() as i32
    };
    let asym = (nl as i32 - nr as i32).unsigned_abs() as i32;
    e += (NINIO_M * asym).min(NINIO_MAX);
    if !(nl == 1 && nr == 1) {
        let b5o = seq[i + 1] as usize;
        let b3o = seq[j - 1] as usize;
        e += INTERIOR_MM[pi][b5o][b3o];
        let b3i = if p > 0 { seq[p - 1] as usize } else { 0 };
        let b5i = if q + 1 < seq.len() { seq[q + 1] as usize } else { 0 };
        e += INTERIOR_MM[pp][b3i][b5i];
    }
    if is_au_gu(pi) { e += TERMINAL_AU_PENALTY; }
    if is_au_gu(pp) { e += TERMINAL_AU_PENALTY; }
    e as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tetraloop_energy_no_double_mismatch() {
        // GCCUCCGGGC with structure ((........))
        // CUCCGG is a special tetraloop with bonus 270
        // Encoding: G=2, C=1, U=3
        let _seq: Vec<u8> = vec![2, 1, 1, 3, 1, 1, 2, 2, 2, 1];
        //                      G  C  C  U  C  C  G  G  G  C
        //                      0  1  2  3  4  5  6  7  8  9
        
        // Hairpin from position 1 to 8 (C-G closing pair)
        // Loop content: CUCCGG (positions 1..=8) → size = 8-1-1 = 6
        // That's NOT a tetraloop (size != 4).
        //
        // For a real tetraloop test, we need size == 4.
        // CCUCCGG: C(1) closing with G(8), loop = UCCG, size = 4
        // The window seq[1..=6] = [1,3,1,1,2,2] = CUCCGG → bonus 270
        //
        // Use GCUCCGGC: G(0)-C(7), inner C(1)-G(6), hairpin at (1,6), size=4
        // seq[1..=6] = [1,3,1,1,2,2] = CUCCGG
        let seq2: Vec<u8> = vec![2, 1, 3, 1, 1, 2, 2, 1];
        //                       G  C  U  C  C  G  G  C
        //                       0  1  2  3  4  5  6  7
        // Hairpin at (1, 6): closing pair C(1)-G(6), size = 6-1-1 = 4
        // Window seq[1..=6] = [1,3,1,1,2,2] = CUCCGG → bonus 270
        let pi = pair_index(seq2[1], seq2[6]).unwrap(); // C-G = 2
        let energy = hairpin_e(&seq2, 1, 6, pi);
        
        // Should be exactly 270 (the tetraloop bonus)
        // NOT 270 + mismatch
        assert_eq!(energy, 270, 
            "Tetraloop energy should be exactly 270, got {}", energy);
    }
    
    #[test]
    fn test_dp_matches_eval_energy() {
        // After fixing eval_combined.rs, verify consistency
        let seq: Vec<u8> = vec![2, 2, 2, 2, 0, 0, 0, 0, 1, 1, 1, 1];
        //                      G  G  G  G  A  A  A  A  C  C  C  C
        
        let (dp_energy, pairs, _structure) = fold_mfe_full(&seq);
        let eval_e = crate::energy::eval_energy(&seq, &pairs);
        
        assert_eq!(dp_energy, eval_e,
            "DP energy ({}) should match eval_energy ({})", 
            dp_energy, eval_e);
    }

    #[test]
    fn test_gaaa_tetraloop_folds_correctly() {
        // GGGGAAAACCCC should fold to ((((....))))
        let seq: Vec<u8> = vec![2, 2, 2, 2, 0, 0, 0, 0, 1, 1, 1, 1];
        let (mfe, _pairs, structure) = fold_mfe_full(&seq);
        
        assert_eq!(structure, "((((....))))");
        // Energy should be reasonable (not double-counted)
        assert!(mfe > -600, "MFE {} suspiciously low", mfe);
    }

    #[test]
    fn test_dp_eval_parity_multiple_sequences() {
        let sequences: Vec<Vec<u8>> = vec![
            vec![2, 2, 2, 0, 0, 0, 1, 1, 1],           // GGGAAACCC
            vec![2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1],   // GCGCGCGCGCGC
            vec![0, 3, 1, 2, 0, 3, 1, 2, 0, 3, 1, 2, 0, 3, 1, 2], // AUCGAUCGAUCGAUCG
            vec![1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2], // CCCCGGGGAAAACCCCGGGG
        ];
        
        for seq in &sequences {
            let (dp_e, pairs, db) = fold_mfe_full(seq);
            let eval_e = crate::energy::eval_energy(seq, &pairs);
            assert_eq!(dp_e, eval_e,
                "DP ({}) != eval ({}) for structure {}", dp_e, eval_e, db);
        }
    }
}
