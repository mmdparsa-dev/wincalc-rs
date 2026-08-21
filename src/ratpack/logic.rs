use crate::ratpack::*;
use crate::ratpack::types::*;
use std::cmp::{max, min};

pub fn lshrat(pa: &mut PRat, b: &PRat, radix: u32, precision: i32) {
    intrat(pa, radix, precision);
    if !zernum(&pa.pp) {
        if rat_gt(b, &rat_max_exp(), precision) {
            panic!("CALC_E_DOMAIN");
        }
        let intb = rattoi32(b, radix, precision);
        let mut pwr = duprat(&rat_two());
        ratpowi32(&mut pwr, intb, precision);
        mulrat(pa, &pwr, precision);
    }
}

pub fn rshrat(pa: &mut PRat, b: &PRat, radix: u32, precision: i32) {
    intrat(pa, radix, precision);
    if !zernum(&pa.pp) {
        if rat_lt(b, &rat_min_exp(), precision) {
            panic!("CALC_E_DOMAIN");
        }
        let intb = rattoi32(b, radix, precision);
        let mut pwr = duprat(&rat_two());
        ratpowi32(&mut pwr, intb, precision);
        divrat(pa, &pwr, precision);
    }
}

#[derive(Clone, Copy)]
pub enum BoolFuncs {
    FuncAnd,
    FuncOr,
    FuncXor,
}

pub fn andrat(pa: &mut PRat, b: &PRat, radix: u32, precision: i32) {
    boolrat(pa, b, BoolFuncs::FuncAnd, radix, precision);
}

pub fn orrat(pa: &mut PRat, b: &PRat, radix: u32, precision: i32) {
    boolrat(pa, b, BoolFuncs::FuncOr, radix, precision);
}

pub fn xorrat(pa: &mut PRat, b: &PRat, radix: u32, precision: i32) {
    boolrat(pa, b, BoolFuncs::FuncXor, radix, precision);
}

fn boolrat(pa: &mut PRat, b: &PRat, func: BoolFuncs, radix: u32, precision: i32) {
    intrat(pa, radix, precision);
    let mut tmp = duprat(b);
    intrat(&mut tmp, radix, precision);

    boolnum(&mut pa.pp, &tmp.pp, func);
}

fn boolnum(pa: &mut PNumber, b: &PNumber, func: BoolFuncs) {
    let a = pa.clone();
    let mut cdigits = max(a.cdigit + a.exp, b.cdigit + b.exp) - min(a.exp, b.exp);
    let mut c = createnum(cdigits as u32);
    c.exp = min(a.exp, b.exp);
    let mut mexp = c.exp;
    c.cdigit = cdigits;
    
    let orig_c_cdigit = c.cdigit;
    
    let mut pcha_idx = 0;
    let mut pchb_idx = 0;
    let mut pchc_idx = 0;
    
    while cdigits > 0 {
        let da = if mexp >= a.exp && (cdigits + a.exp - c.exp > (orig_c_cdigit - a.cdigit)) {
            let val = a.mant.get(pcha_idx).copied().unwrap_or(0);
            pcha_idx += 1;
            val
        } else {
            0
        };
        
        let db = if mexp >= b.exp && (cdigits + b.exp - c.exp > (orig_c_cdigit - b.cdigit)) {
            let val = b.mant.get(pchb_idx).copied().unwrap_or(0);
            pchb_idx += 1;
            val
        } else {
            0
        };
        
        let res = match func {
            BoolFuncs::FuncAnd => da & db,
            BoolFuncs::FuncOr => da | db,
            BoolFuncs::FuncXor => da ^ db,
        };
        
        if pchc_idx < c.mant.len() {
            c.mant[pchc_idx] = res;
        } else {
            c.mant.push(res);
        }
        pchc_idx += 1;
        
        cdigits -= 1;
        mexp += 1;
    }
    
    c.sign = a.sign;
    
    let mut i = c.cdigit as usize;
    while i > 1 && c.mant.get(i - 1).copied().unwrap_or(1) == 0 {
        c.cdigit -= 1;
        i -= 1;
    }
    
    *pa = c;
}

pub fn remrat(pa: &mut PRat, b: &PRat) {
    if zerrat(b) {
        panic!("CALC_E_INDEFINITE");
    }

    let mut tmp = duprat(b);

    mulnumx(&mut pa.pp, &tmp.pq);
    mulnumx(&mut tmp.pp, &pa.pq);
    remnum(&mut pa.pp, &tmp.pp, BASEX);
    mulnumx(&mut pa.pq, &tmp.pq);

    renormalize(pa);
}

pub fn modrat(pa: &mut PRat, b: &PRat) {
    if zerrat(b) {
        return;
    }

    let mut tmp = duprat(b);

    let sign_pa = pa.pp.sign * pa.pq.sign;
    let sign_b = b.pp.sign * b.pq.sign;
    
    let need_adjust = if sign_pa == -1 {
        sign_b == 1
    } else {
        sign_b == -1
    };

    mulnumx(&mut pa.pp, &tmp.pq);
    mulnumx(&mut tmp.pp, &pa.pq);
    remnum(&mut pa.pp, &tmp.pp, BASEX);
    mulnumx(&mut pa.pq, &tmp.pq);

    if need_adjust && !zerrat(pa) {
        _addrat(pa, b, BASEX);
    }

    renormalize(pa);
}
