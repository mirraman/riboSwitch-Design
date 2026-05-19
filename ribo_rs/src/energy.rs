use crate::params::*;
pub fn eval_energy(seq: &[u8], pair_table: &[i32]) -> i32 {
    let n = seq.len();
    assert_eq!(n, pair_table.len());
    let mut total = 0i32;
    let mut top_pairs: Vec<(usize, usize)> = Vec::new();
    let mut pos = 0usize;
    while pos < n {
        let j = pair_table[pos];
        if j > pos as i32 {
            top_pairs.push((pos, j as usize));
            pos = j as usize + 1;
        } else {
            pos += 1;
        }
    }
    total += external_energy(seq, pair_table, &top_pairs, n);
    for &(i, j) in &top_pairs {
        total += process_pair(seq, pair_table, i, j);
    }
    total
}
fn process_pair(seq: &[u8], pair_table: &[i32], i: usize, j: usize) -> i32 {
    let children = find_children(pair_table, i, j);
    let mut total = match children.len() {
        0 => hairpin_energy(seq, i, j),
        1 => {
            let (p, q) = children[0];
            let nl = p - i - 1;
            let nr = j - q - 1;
            if nl == 0 && nr == 0 {
                stack_energy(seq, i, j, p, q)
            } else if nl == 0 || nr == 0 {
                bulge_energy(seq, i, j, p, q, nl + nr)
            } else {
                interior_energy(seq, i, j, p, q, nl, nr)
            }
        }
        _ => multiloop_energy(seq, pair_table, i, j, &children),
    };
    for &(p, q) in &children {
        total += process_pair(seq, pair_table, p, q);
    }
    total
}
fn find_children(pair_table: &[i32], i: usize, j: usize) -> Vec<(usize, usize)> {
    let mut children = Vec::new();
    let mut pos = i + 1;
    while pos < j {
        let q = pair_table[pos];
        if q > pos as i32 && (q as usize) < j {
            children.push((pos, q as usize));
            pos = q as usize + 1;
        } else {
            pos += 1;
        }
    }
    children
}
fn hairpin_energy(seq: &[u8], i: usize, j: usize) -> i32 {
    let size = j - i - 1;
    if size < MIN_HAIRPIN {
        return INF;
    }
    let pi = match pair_index(seq[i], seq[j]) {
        Some(p) => p,
        None => return INF,
    };
    let mut energy = if size <= 30 {
        HAIRPIN_INIT[size]
    } else {
        HAIRPIN_INIT[30]
            + (LOOP_EXTRAPOLATION_COEFF * (size as f64 / 30.0).ln() * 100.0).round() as i32
    };
    if size == 3 {
        let w = &seq[i..=j];
        let bonus = triloop_bonus(w);
        if bonus != 0 {
            energy = bonus;
        }
        if is_au_gu(pi) {
            energy += TERMINAL_AU_PENALTY;
        }
        if (i + 1..j).all(|k| seq[k] == 1) {
            energy += HAIRPIN_C3;
        }
        return energy;
    }
    if size == 4 {
        let w = &seq[i..=j];
        let bonus = tetraloop_bonus(w);
        if bonus != 0 {
            return bonus;
        }
    }
    let b5 = seq[i + 1] as usize;
    let b3 = seq[j - 1] as usize;
    energy += HAIRPIN_MM[pi][b5][b3];
    let mm = (seq[i + 1], seq[j - 1]);
    if mm == (3, 3) || mm == (2, 0) {
        energy += HAIRPIN_UU_GA_BONUS;
    }
    if mm == (2, 2) {
        energy += HAIRPIN_GG_BONUS;
    }
    if pi == 4 && i >= 2 && seq[i - 1] == 2 && seq[i - 2] == 2 {
        energy += HAIRPIN_SPECIAL_GU;
    }
    if (i + 1..j).all(|k| seq[k] == 1) {
        energy += HAIRPIN_C_SLOPE * size as i32 + HAIRPIN_C_INTERCEPT;
    }
    energy
}
fn stack_energy(seq: &[u8], i: usize, j: usize, p: usize, q: usize) -> i32 {
    let oi = match pair_index(seq[i], seq[j]) { Some(x) => x, None => return INF };
    let ii = match pair_index(seq[p], seq[q]) { Some(x) => x, None => return INF };
    STACK[oi][ii]
}
fn bulge_energy(seq: &[u8], i: usize, j: usize, p: usize, q: usize, size: usize) -> i32 {
    let mut energy = if size <= 30 {
        BULGE_INIT[size]
    } else {
        BULGE_INIT[30]
            + (LOOP_EXTRAPOLATION_COEFF * (size as f64 / 30.0).ln() * 100.0).round() as i32
    };
    let oi = match pair_index(seq[i], seq[j]) { Some(x) => x, None => return INF };
    let ii = match pair_index(seq[p], seq[q]) { Some(x) => x, None => return INF };
    if size == 1 {
        energy += STACK[oi][ii];
    }
    if is_au_gu(oi) { energy += TERMINAL_AU_PENALTY; }
    if is_au_gu(ii) { energy += TERMINAL_AU_PENALTY; }
    energy
}
fn interior_energy(
    seq: &[u8],
    i: usize, j: usize,
    p: usize, q: usize,
    nl: usize, nr: usize,
) -> i32 {
    let oi = match pair_index(seq[i], seq[j]) { Some(x) => x, None => return INF };
    let ii = match pair_index(seq[p], seq[q]) { Some(x) => x, None => return INF };
    let total = nl + nr;
    if nl == 1 && nr == 1 {
        let mm5 = seq[i + 1] as usize;
        let mm3 = seq[j - 1] as usize;
        let val = INT11[oi][ii][mm5][mm3];
        if val < INF {
            return val;
        }
    }
    if nl == 1 && nr == 2 {
        let mm5 = seq[i + 1] as usize;
        let mm3 = seq[j - 1] as usize;
        let mm_mid = seq[j - 2] as usize;
        let val = INT21[oi][ii][mm5][mm3][mm_mid];
        if val < INF {
            return val;
        }
    }
    if nl == 2 && nr == 1 {
        let mm5 = seq[i + 1] as usize;
        let mm3 = seq[j - 1] as usize;
        let mm_mid = seq[i + 2] as usize;
        let val = INT12[oi][ii][mm5][mm3][mm_mid];
        if val < INF {
            return val;
        }
    }
    if nl == 2 && nr == 2 {
        let m5o = seq[i + 1] as usize;
        let m3o = seq[j - 1] as usize;
        let m5i = seq[i + 2] as usize;
        let m3i = seq[j - 2] as usize;
        let val = INT22[oi][ii][m5o][m3o][m5i][m3i];
        if val < INF {
            return val;
        }
    }
    let mut energy = if total <= 30 {
        INTERIOR_INIT[total]
    } else {
        INTERIOR_INIT[30]
            + (LOOP_EXTRAPOLATION_COEFF * (total as f64 / 30.0).ln() * 100.0).round() as i32
    };
    let asym = (nl as i32 - nr as i32).unsigned_abs() as i32;
    energy += (NINIO_M * asym).min(NINIO_MAX);
    if !(nl == 1 && nr == 1) {
        let b5o = seq[i + 1] as usize;
        let b3o = seq[j - 1] as usize;
        energy += INTERIOR_MM[oi][b5o][b3o];
        let b3i = if p > 0 { seq[p - 1] as usize } else { 0 };
        let b5i = if q + 1 < seq.len() { seq[q + 1] as usize } else { 0 };
        energy += INTERIOR_MM[ii][b3i][b5i];
    }
    if is_au_gu(oi) { energy += TERMINAL_AU_PENALTY; }
    if is_au_gu(ii) { energy += TERMINAL_AU_PENALTY; }
    energy
}
fn multiloop_energy(
    seq: &[u8],
    _pair_table: &[i32],
    i: usize, j: usize,
    children: &[(usize, usize)],
) -> i32 {
    let nb = (children.len() as i32) + 1;
    let mut n_unpaired = 0i32;
    let mut prev = i;
    for &(p, q) in children {
        n_unpaired += (p - prev - 1) as i32;
        prev = q;
    }
    n_unpaired += (j - prev - 1) as i32;
    
    // Energia de bază multi-loop
    let mut energy = ML_OFFSET + ML_PER_BRANCH * nb + ML_PER_UNPAIRED * n_unpaired;
    
    // Penalizare AU/GU pentru perechea de închidere
    let ci = match pair_index(seq[i], seq[j]) { Some(x) => x, None => return INF };
    if is_au_gu(ci) { 
        energy += TERMINAL_AU_PENALTY; 
    }
    
    // ========================================================================
    // DANGLES pentru perechea de închidere (looking INTO the multiloop)
    // ========================================================================
    // Dangle 3' de la perechea de închidere (nucleotida i+1)
    if i + 1 < j {
        energy += DANGLE3[ci][seq[i + 1] as usize];
    }
    
    // Dangle 5' de la perechea de închidere (nucleotida j-1)
    if j > i + 1 {
        energy += DANGLE5[ci][seq[j - 1] as usize];
    }
    
    // ========================================================================
    // Procesare ramuri (children)
    // ========================================================================
    for (idx, &(p, q)) in children.iter().enumerate() {
        let bi = match pair_index(seq[p], seq[q]) { Some(x) => x, None => continue };
        
        // Penalizare AU/GU pentru această ramură
        if is_au_gu(bi) { 
            energy += TERMINAL_AU_PENALTY; 
        }
        
        // ====================================================================
        // DANGLES pentru această ramură
        // ====================================================================
        
        // Dangle 5' (nucleotida ÎNAINTE de p)
        if p > i + 1 {
            energy += DANGLE5[bi][seq[p - 1] as usize];
        }
        
        // Dangle 3' (nucleotida DUPĂ q)
        if q + 1 < j {
            if q + 1 < seq.len() {
                energy += DANGLE3[bi][seq[q + 1] as usize];
            }
        }
    }
    
    energy
}
fn external_energy(
    seq: &[u8],
    _pair_table: &[i32],
    top_pairs: &[(usize, usize)],
    n: usize,
) -> i32 {
    let mut energy = 0i32;
    for &(i, j) in top_pairs {
        let pi = match pair_index(seq[i], seq[j]) { Some(x) => x, None => continue };
        if is_au_gu(pi) {
            energy += TERMINAL_AU_PENALTY;
        }
        if i > 0 {
            energy += DANGLE5[pi][seq[i - 1] as usize];
        }
        if j + 1 < n {
            energy += DANGLE3[pi][seq[j + 1] as usize];
        }
    }
    energy
}
