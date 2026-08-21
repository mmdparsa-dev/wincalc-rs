// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `support.cpp`.
//! Contains support functions for rationals and numbers.

use std::cell::RefCell;
use std::rc::Rc;

use crate::ratpack::errors::{CalcError, CalcResult};
use crate::ratpack::types::{AngleType, Number, PNumber, PRat, Rat};

// Global precision tracking (mirrors cbitsofprecision in support.cpp)
thread_local! {
    static CBITSOFPRECISION: RefCell<i32> = RefCell::new(0);
    static G_FTRUEINFINITE: RefCell<bool> = RefCell::new(false);
    /// g_ratio: digits in current radix within BASEX; used for length calculations
    pub static G_RATIO: RefCell<i32> = RefCell::new(0);
}

pub fn get_g_ratio() -> i32 {
    G_RATIO.with(|r| *r.borrow())
}

pub fn get_g_ftrueinfinite() -> bool {
    G_FTRUEINFINITE.with(|r| *r.borrow())
}

pub fn set_g_ftrueinfinite(val: bool) {
    G_FTRUEINFINITE.with(|r| *r.borrow_mut() = val);
}

/// Helper: create a PNumber with a given sign, exp, and mantissa digits.
pub fn make_number(sign: i32, exp: i32, mant: Vec<u32>) -> PNumber {
    Rc::new(RefCell::new(Number { sign, exp, mant }))
}

/// Helper: create a PRat from two PNumbers.
pub fn make_rat(pp: PNumber, pq: PNumber) -> PRat {
    Rc::new(RefCell::new(Rat { pp, pq }))
}

/// dupnum: deep-clone a Number, replacing the destination.
pub fn dupnum(dest: &mut PNumber, src: &PNumber) {
    let cloned = src.borrow().clone();
    *dest = Rc::new(RefCell::new(cloned));
}

/// duprat: deep-clone a Rat, replacing the destination.
pub fn duprat(dest: &mut PRat, src: &PRat) {
    let src_b = src.borrow();
    let new_pp = {
        let n = src_b.pp.borrow().clone();
        Rc::new(RefCell::new(n))
    };
    let new_pq = {
        let n = src_b.pq.borrow().clone();
        Rc::new(RefCell::new(n))
    };
    *dest = Rc::new(RefCell::new(Rat {
        pp: new_pp,
        pq: new_pq,
    }));
}

/// zernum: returns true if the number's first mantissa digit is 0 (or mant is empty).
pub fn zernum(n: &PNumber) -> bool {
    let b = n.borrow();
    b.mant.is_empty() || b.mant[0] == 0
}

/// zerrat: returns true if the numerator is zero.
pub fn zerrat(r: &PRat) -> bool {
    zernum(&r.borrow().pp)
}

/// equnum: returns true if two Numbers are equal in value.
pub fn equnum(a: &PNumber, b: &PNumber) -> bool {
    let ab = a.borrow();
    let bb = b.borrow();
    ab.sign == bb.sign && ab.exp == bb.exp && ab.mant == bb.mant
}

/// sign of rational: sign(pp) * sign(pq)
pub fn rat_sign(r: &PRat) -> i32 {
    let b = r.borrow();
    b.pp.borrow().sign * b.pq.borrow().sign
}

/// absrat: sets both pp and pq signs to +1 (absolute value of rational in-place).
pub fn absrat(r: &PRat) {
    let b = r.borrow();
    b.pp.borrow_mut().sign = 1;
    b.pq.borrow_mut().sign = 1;
}

/// lognum2: integral portion of log base-BASEX of a number.
pub fn lognum2(n: &PNumber) -> i32 {
    let b = n.borrow();
    b.mant.len() as i32 + b.exp
}

/// lograt2: difference of lognum2 for pp and pq.
pub fn lograt2(r: &PRat) -> i32 {
    let b = r.borrow();
    lognum2(&b.pp) - lognum2(&b.pq)
}

/// lognumradix: approximate log in current radix
pub fn lognumradix(n: &PNumber) -> i32 {
    let b = n.borrow();
    (b.mant.len() as i32 + b.exp) * get_g_ratio()
}

/// logratradix: difference of lognumradix for pp and pq
pub fn logratradix(r: &PRat) -> i32 {
    let b = r.borrow();
    lognumradix(&b.pp) - lognumradix(&b.pq)
}

// ---------------------------------------------------------------------------
//   Comparison functions (rat_equ, rat_gt, rat_lt, rat_ge, rat_le, rat_neq)
//   These all work by subtracting b from a (using _addrat with negated sign)
//   and testing the sign of the result.
// ---------------------------------------------------------------------------

pub fn rat_equ(a: &PRat, b: &PRat, precision: i32) -> CalcResult<bool> {
    let mut rattmp = Rc::clone(a);
    // Temporarily negate numerator of a (clone first)
    duprat(&mut rattmp, a);
    rattmp.borrow().pp.borrow_mut().sign *= -1;
    _addrat(&mut rattmp, b, precision)?;
    let bret = zernum(&rattmp.borrow().pp);
    Ok(bret)
}

pub fn rat_neq(a: &PRat, b: &PRat, precision: i32) -> CalcResult<bool> {
    Ok(!rat_equ(a, b, precision)?)
}

pub fn rat_gt(a: &PRat, b: &PRat, precision: i32) -> CalcResult<bool> {
    let mut rattmp = Rc::clone(a);
    duprat(&mut rattmp, a);
    // Temporarily negate sign of b's numerator (without mutating b permanently)
    b.borrow().pp.borrow_mut().sign *= -1;
    _addrat(&mut rattmp, b, precision)?;
    b.borrow().pp.borrow_mut().sign *= -1;
    let bret = !zernum(&rattmp.borrow().pp) && rat_sign(&rattmp) == 1;
    Ok(bret)
}

pub fn rat_ge(a: &PRat, b: &PRat, precision: i32) -> CalcResult<bool> {
    let mut rattmp = Rc::clone(a);
    duprat(&mut rattmp, a);
    b.borrow().pp.borrow_mut().sign *= -1;
    _addrat(&mut rattmp, b, precision)?;
    b.borrow().pp.borrow_mut().sign *= -1;
    let bret = zernum(&rattmp.borrow().pp) || rat_sign(&rattmp) == 1;
    Ok(bret)
}

pub fn rat_lt(a: &PRat, b: &PRat, precision: i32) -> CalcResult<bool> {
    let mut rattmp = Rc::clone(a);
    duprat(&mut rattmp, a);
    b.borrow().pp.borrow_mut().sign *= -1;
    _addrat(&mut rattmp, b, precision)?;
    b.borrow().pp.borrow_mut().sign *= -1;
    let bret = !zernum(&rattmp.borrow().pp) && rat_sign(&rattmp) == -1;
    Ok(bret)
}

pub fn rat_le(a: &PRat, b: &PRat, precision: i32) -> CalcResult<bool> {
    let mut rattmp = Rc::clone(a);
    duprat(&mut rattmp, a);
    b.borrow().pp.borrow_mut().sign *= -1;
    _addrat(&mut rattmp, b, precision)?;
    b.borrow().pp.borrow_mut().sign *= -1;
    let bret = zernum(&rattmp.borrow().pp) || rat_sign(&rattmp) == -1;
    Ok(bret)
}

// ---------------------------------------------------------------------------
//   Forward declarations (implemented in other modules of crate::ratpack)
// ---------------------------------------------------------------------------
//
// The following functions are declared here but implemented in their respective
// modules (rat.rs, num.rs, conv.rs, etc.). They are brought in via `pub use`.
// We use function pointers/stubs here to break potential circular dependencies.
// The real implementations are linked at compile time via the module tree.

/// _addrat stub: implemented in rat.rs
pub fn _addrat(pa: &mut PRat, b: &PRat, precision: i32) -> CalcResult<()> {
    // Calls into rat module -- linked by module inclusion in mod.rs
    crate::ratpack::rat::_addrat(pa, b, precision)
}

// ---------------------------------------------------------------------------
//   intrat: truncates a rational to its integer part.
// ---------------------------------------------------------------------------

pub fn intrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    // Only truncate if nonzero and denominator isn't 1.
    let is_zero = zernum(&px.borrow().pp);
    let denom_is_one = {
        let num_one = crate::ratpack::constants::num_one();
        equnum(&px.borrow().pq, &num_one)
    };

    if !is_zero && !denom_is_one {
        crate::ratpack::conv::flatrat(px, radix, precision)?;

        let mut pret = Rc::clone(px);
        duprat(&mut pret, px);
        crate::ratpack::rat::remrat(
            &mut pret,
            &crate::ratpack::constants::rat_one(),
        )?;

        // Flatten pret if denominators differ after remrat
        let denoms_differ = {
            let pb = px.borrow();
            let rb = pret.borrow();
            !equnum(&pb.pq, &rb.pq)
        };
        if denoms_differ {
            crate::ratpack::conv::flatrat(&mut pret, radix, precision)?;
        }

        crate::ratpack::rat::_subrat(px, &pret, precision)?;
        crate::ratpack::conv::flatrat(px, radix, precision)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
//   scale: scales x to within the range [0, scalefact)
// ---------------------------------------------------------------------------

pub fn scale(px: &mut PRat, scalefact: &PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let mut pret: PRat = Rc::clone(px);
    duprat(&mut pret, px);

    // logscale: extra precision needed
    let logscale = {
        let b = pret.borrow();
        let pp = b.pp.borrow();
        let pq = b.pq.borrow();
        get_g_ratio()
            * ((pp.mant.len() as i32 + pp.exp) - (pq.mant.len() as i32 + pq.exp))
    };
    let prec = if logscale > 0 {
        precision + logscale
    } else {
        precision
    };

    crate::ratpack::rat::divrat(&mut pret, scalefact, prec)?;
    intrat(&mut pret, radix, prec)?;
    crate::ratpack::rat::mulrat(&mut pret, scalefact, prec)?;
    pret.borrow().pp.borrow_mut().sign *= -1;
    _addrat(px, &pret, prec)?;
    Ok(())
}

// ---------------------------------------------------------------------------
//   scale2pi: scales x to [0, 2π)
// ---------------------------------------------------------------------------

pub fn scale2pi(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    use crate::ratpack::constants::{rat_half, rat_six, rat_two};
    use crate::ratpack::itrans::asinrat;
    use crate::ratpack::rat::{divrat, mulrat};

    let mut pret: PRat = Rc::clone(px);
    duprat(&mut pret, px);

    let logscale = {
        let b = pret.borrow();
        let pp = b.pp.borrow();
        let pq = b.pq.borrow();
        get_g_ratio()
            * ((pp.mant.len() as i32 + pp.exp) - (pq.mant.len() as i32 + pq.exp))
    };

    let (prec, mut my_two_pi) = if logscale > 0 {
        let new_prec = precision + logscale;
        let mut tpi = rat_half();
        asinrat(&mut tpi, radix, new_prec)?;
        mulrat(&mut tpi, &rat_six(), new_prec)?;
        mulrat(&mut tpi, &rat_two(), new_prec)?;
        (new_prec, tpi)
    } else {
        (precision, crate::ratpack::support::TWO_PI.with(|v| Rc::clone(&v.borrow())))
    };

    divrat(&mut pret, &my_two_pi, prec)?;
    intrat(&mut pret, radix, prec)?;
    mulrat(&mut pret, &my_two_pi, prec)?;
    pret.borrow().pp.borrow_mut().sign *= -1;
    _addrat(px, &pret, prec)?;
    Ok(())
}

// ---------------------------------------------------------------------------
//   inbetween: clamps *px to [-range, +range]
// ---------------------------------------------------------------------------

pub fn inbetween(px: &mut PRat, range: &PRat, precision: i32) -> CalcResult<()> {
    if rat_gt(px, range, precision)? {
        duprat(px, range);
    } else {
        range.borrow().pp.borrow_mut().sign *= -1;
        if rat_lt(px, range, precision)? {
            duprat(px, range);
        }
        range.borrow().pp.borrow_mut().sign *= -1;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
//   trimit: chops digits from rational to control time complexity.
// ---------------------------------------------------------------------------

pub fn trimit(px: &mut PRat, precision: i32) -> CalcResult<()> {
    if get_g_ftrueinfinite() {
        return Ok(());
    }

    let trim_val = {
        let b = px.borrow();
        let pp = b.pp.borrow();
        let pq = b.pq.borrow();
        let pp_size = pp.mant.len() as i32 + pp.exp;
        let pq_size = pq.mant.len() as i32 + pq.exp;
        get_g_ratio() * (pp_size.min(pq_size) - 1) - precision
    };

    if trim_val > get_g_ratio() {
        let trim = trim_val / get_g_ratio();

        // Trim pp
        {
            let b = px.borrow();
            let mut pp = b.pp.borrow_mut();
            if trim <= pp.exp {
                pp.exp -= trim;
            } else {
                let offset = (trim - pp.exp) as usize;
                let new_len = pp.mant.len().saturating_sub(offset);
                pp.mant.drain(..offset.min(pp.mant.len()));
                pp.mant.truncate(new_len);
                pp.exp = 0;
            }
        }

        // Trim pq
        {
            let b = px.borrow();
            let mut pq = b.pq.borrow_mut();
            if trim <= pq.exp {
                pq.exp -= trim;
            } else {
                let offset = (trim - pq.exp) as usize;
                let new_len = pq.mant.len().saturating_sub(offset);
                pq.mant.drain(..offset.min(pq.mant.len()));
                pq.mant.truncate(new_len);
                pq.exp = 0;
            }
        }
    }

    // Normalize exponents
    let min_exp = {
        let b = px.borrow();
        b.pp.borrow().exp.min(b.pq.borrow().exp)
    };
    px.borrow().pp.borrow_mut().exp -= min_exp;
    px.borrow().pq.borrow_mut().exp -= min_exp;

    Ok(())
}

// ---------------------------------------------------------------------------
//   Thread-local globals for mathematical constants.
//   These are initialized by change_constants() and consumed by the engine.
// ---------------------------------------------------------------------------
//
//   In C++, these were mutable globals initialized lazily in ChangeConstants().
//   In Rust we use thread_local! with RefCell to get the same semantics safely.

use std::cell::RefCell as TLRefCell;

thread_local! {
    pub static TWO_PI: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static PI_OVER_TWO: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static ONE_PT_FIVE_PI: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static E_TO_ONE_HALF: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static RAT_EXP: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static RAD_TO_DEG: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static RAD_TO_GRAD: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static RAT_NRADIX: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static RAT_SMALLEST: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
    pub static RAT_NEGSMALLEST: TLRefCell<PRat> = TLRefCell::new(crate::ratpack::constants::rat_zero());
}
