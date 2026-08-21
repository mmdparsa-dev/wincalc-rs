// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `exp.cpp`.
//! Contains exp, log, log10, pow, and related functions for rationals.

use crate::ratpack::constants::{
    rat_half, rat_max_exp, rat_min_exp, rat_negsmallest, rat_one, rat_smallest, rat_ten, rat_two,
    rat_zero,
};
use crate::ratpack::conv::i32tonum;
use crate::ratpack::errors::{CalcError, CalcResult};
use crate::ratpack::rat::{_addrat, _subrat, divrat, fracrat, mulrat};
use crate::ratpack::support::{
    absrat, dupnum, duprat, equnum, get_g_ratio, inbetween, intrat, lograt2, rat_gt, rat_lt,
    rat_sign, trimit, zernum, zerrat,
};
use crate::ratpack::types::PRat;
use std::rc::Rc;

// ---------------------------------------------------------------------------
//  Taylor series helpers
// ---------------------------------------------------------------------------
//
//  CREATETAYLOR / DESTROYTAYLOR / NEXTTERM / SMALL_ENOUGH_RAT were macros in C++.
//  We represent the Taylor state as a local struct to avoid any unsafe sharing.

struct TaylorState {
    /// Return accumulator (pret in C++)
    pret: PRat,
    /// Current term (thisterm in C++)
    thisterm: PRat,
    /// x^2 factor driving the series (xx in C++)
    xx: PRat,
    /// Counter for term index (n2 in C++)
    n2: crate::ratpack::types::PNumber,
}

impl TaylorState {
    /// Creates the Taylor state using the initial x value.
    fn new(x: &PRat) -> Self {
        let mut pret = rat_zero();
        duprat(&mut pret, x);
        let mut thisterm = rat_zero();
        duprat(&mut thisterm, x);
        let mut xx = rat_zero();
        duprat(&mut xx, x);
        let n2 = i32tonum(0, crate::ratpack::types::BASEX);
        TaylorState {
            pret,
            thisterm,
            xx,
            n2,
        }
    }

    /// Returns true when thisterm is small enough relative to pret.
    fn small_enough(&self, precision: i32) -> CalcResult<bool> {
        // SMALL_ENOUGH_RAT: |thisterm| < |pret| * radix^(-precision)
        // Approximated by checking LOGRAT2
        let log_term = crate::ratpack::support::lograt2(&self.thisterm);
        let log_pret = crate::ratpack::support::lograt2(&self.pret);
        Ok(log_term < log_pret - precision)
    }
}

// ---------------------------------------------------------------------------
//  _exprat: core Taylor series for exp(x), for |x| < 1
// ---------------------------------------------------------------------------

pub fn _exprat(px: &mut PRat, precision: i32) -> CalcResult<()> {
    use crate::ratpack::conv::i32tonum;
    use crate::ratpack::num::addnum;
    use crate::ratpack::types::BASEX;

    // pret = 1 (the constant term of exp series)
    let mut pret = rat_zero();
    {
        let one_num = i32tonum(1, BASEX);
        dupnum(&mut pret.borrow_mut().pp, &one_num);
        let one_num2 = i32tonum(1, BASEX);
        dupnum(&mut pret.borrow_mut().pq, &one_num2);
    }
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, &pret);
    let mut xx = rat_zero();
    duprat(&mut xx, px);

    let mut n2 = i32tonum(0, BASEX);

    loop {
        // n2 += 1
        let one = i32tonum(1, BASEX);
        addnum(&mut n2, &one, BASEX);

        // thisterm *= x / n2
        {
            crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &xx.borrow().pp);
            crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &n2);
        }
        trimit(&mut thisterm, precision)?;

        // pret += thisterm
        _addrat(&mut pret, &thisterm, precision)?;

        if small_enough_rat(&thisterm, &pret, precision) {
            break;
        }
    }

    duprat(px, &pret);
    Ok(())
}

fn small_enough_rat(thisterm: &PRat, _pret: &PRat, precision: i32) -> bool {
    crate::ratpack::support::lograt2(thisterm) < -precision
}

// ---------------------------------------------------------------------------
//  exprat: full exp with integer scaling
// ---------------------------------------------------------------------------

pub fn exprat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    use crate::ratpack::conv::{i32torat, rattoi32};
    use crate::ratpack::num::ratpowi32;

    if rat_gt(px, &rat_max_exp(), precision)?
        || rat_lt(px, &rat_min_exp(), precision)?
    {
        return Err(CalcError::Domain);
    }

    // pwr = e (rat_exp constant)
    let mut pwr = crate::ratpack::support::RAT_EXP
        .with(|v| Rc::clone(&v.borrow()));
    let mut e_clone = rat_zero();
    duprat(&mut e_clone, &pwr);

    // pint = floor(x)
    let mut pint = rat_zero();
    duprat(&mut pint, px);
    intrat(&mut pint, radix, precision)?;

    let intpwr = rattoi32(&pint, radix, precision);
    ratpowi32(&mut e_clone, intpwr, precision)?;

    _subrat(px, &pint, precision)?;

    // Check if x is an exact integral power of e (fractional part ≈ 0)
    let frac_is_zero = rat_gt(px, &rat_negsmallest(), precision)?
        && rat_lt(px, &rat_smallest(), precision)?;

    if frac_is_zero {
        duprat(px, &e_clone);
    } else {
        _exprat(px, precision)?;
        mulrat(px, &e_clone, precision)?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
//  __lograt: inner Taylor series for log(x), for x near 1
// ---------------------------------------------------------------------------

fn __lograt(px: &mut PRat, precision: i32) -> CalcResult<()> {
    use crate::ratpack::conv::i32tonum;
    use crate::ratpack::num::addnum;
    use crate::ratpack::types::BASEX;

    // Compute (x-1) first: pq sign flip trick
    px.borrow().pq.borrow_mut().sign *= -1;
    {
        let pq_clone = Rc::clone(&px.borrow().pq);
        let pq_num = pq_clone.borrow().clone();
        let pq_num_rc = Rc::new(std::cell::RefCell::new(pq_num));
        addnum(&mut px.borrow_mut().pp, &pq_num_rc, BASEX);
    }
    px.borrow().pq.borrow_mut().sign *= -1;

    let mut pret = rat_zero();
    duprat(&mut pret, px);
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, px);
    let mut n2 = i32tonum(1, BASEX);
    px.borrow().pp.borrow_mut().sign *= -1;

    loop {
        // thisterm *= x_minus_1 * n2 / (n2+1)
        {
            crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &n2);
            let one = i32tonum(1, BASEX);
            addnum(&mut n2, &one, BASEX);
            crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &n2);
            // multiply by (x-1) direction (px.pp is already negated for series)
            crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &px.borrow().pp);
            crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &px.borrow().pq);
        }
        trimit(&mut thisterm, precision)?;
        let trimmed_log2 = crate::ratpack::support::lograt2(&thisterm);
        let trim_top = trimmed_log2 - get_g_ratio() - precision;
        // TRIMTOP equivalent -- just prune via trimit
        trimit(&mut thisterm, precision)?;

        _addrat(&mut pret, &thisterm, precision)?;

        if small_enough_rat(&thisterm, &pret, precision) {
            break;
        }
    }

    duprat(px, &pret);
    Ok(())
}

// ---------------------------------------------------------------------------
//  _lograt: natural log with scaling
// ---------------------------------------------------------------------------

pub fn _lograt(px: &mut PRat, precision: i32) -> CalcResult<()> {
    use crate::ratpack::conv::i32torat;

    if rat_le(px, &rat_zero(), precision)? {
        return Err(CalcError::Domain);
    }

    let fneglog = rat_lt(px, &rat_one(), precision)?;
    if fneglog {
        // Swap pp and pq to bring x above 1
        let old_pp = Rc::clone(&px.borrow().pp);
        let old_pq = Rc::clone(&px.borrow().pq);
        px.borrow_mut().pp = old_pq;
        px.borrow_mut().pq = old_pp;
    }

    let mut pwr = rat_zero();
    // Scale by BASEXPWR powers of 2
    let log2_x = lograt2(px);
    if log2_x > 1 {
        let intpwr = log2_x - 1;
        px.borrow().pq.borrow_mut().exp += intpwr;
        let mut pwr2 = i32torat(intpwr * crate::ratpack::types::BASEXPWR as i32);
        mulrat(&mut pwr2, &crate::ratpack::support::LN_TWO.with(|v| Rc::clone(&v.borrow())), precision)?;
        duprat(&mut pwr, &pwr2);
        trimit(px, precision)?;
    } else {
        duprat(&mut pwr, &rat_zero());
    }

    let mut offset = rat_zero();
    duprat(&mut offset, &rat_zero());

    let e_half = crate::ratpack::support::E_TO_ONE_HALF.with(|v| Rc::clone(&v.borrow()));
    while rat_gt(px, &e_half, precision)? {
        divrat(px, &e_half, precision)?;
        _addrat(&mut offset, &rat_one(), precision)?;
    }

    __lograt(px, precision)?;

    // offset /= 2; pwr += offset
    divrat(&mut offset, &rat_two(), precision)?;
    _addrat(&mut pwr, &offset, precision)?;
    _addrat(px, &pwr, precision)?;

    trimit(px, precision)?;

    if fneglog {
        px.borrow().pp.borrow_mut().sign *= -1;
    }

    Ok(())
}

pub fn lograt(px: &mut PRat, precision: i32) -> CalcResult<()> {
    use crate::ratpack::rat::_snaprat;
    let mut a = rat_zero();
    duprat(&mut a, px);
    _lograt(px, precision)?;
    _snaprat(px, &a, None, precision)?;
    Ok(())
}

pub fn log10rat(px: &mut PRat, precision: i32) -> CalcResult<()> {
    lograt(px, precision)?;
    let ln_ten = crate::ratpack::support::LN_TEN.with(|v| Rc::clone(&v.borrow()));
    divrat(px, &ln_ten, precision)?;
    Ok(())
}

// ---------------------------------------------------------------------------
//  IsEven: returns true if numerator of x (with denom 1) is even
// ---------------------------------------------------------------------------

pub fn is_even(x: &PRat, radix: u32, precision: i32) -> CalcResult<bool> {
    let mut tmp = rat_zero();
    duprat(&mut tmp, x);
    divrat(&mut tmp, &rat_two(), precision)?;
    fracrat(&mut tmp, radix, precision)?;
    _addrat(&mut tmp, &tmp.clone(), precision)?;
    _subrat(&mut tmp, &rat_one(), precision)?;
    Ok(rat_lt(&tmp, &rat_zero(), precision)?)
}

// ---------------------------------------------------------------------------
//  ratpowi32: x^n for integer n (from num.rs, referenced here)
//  rootrat: x^(1/y) -- ln(x)*y then exp
// ---------------------------------------------------------------------------

pub fn rootrat(px: &mut PRat, y: &PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let inv_y = {
        let mut r = rat_one();
        divrat(&mut r, y, precision)?;
        r
    };
    powratcomp(px, &inv_y, radix, precision)
}

// ---------------------------------------------------------------------------
//  powratcomp: core x^y implementation
// ---------------------------------------------------------------------------

pub fn powratcomp(px: &mut PRat, y: &PRat, radix: u32, precision: i32) -> CalcResult<()> {
    use crate::ratpack::conv::rattoi32;
    use crate::ratpack::num::ratpowi32;

    let sign = rat_sign(px);
    px.borrow().pp.borrow_mut().sign = 1;
    px.borrow().pq.borrow_mut().sign = 1;

    if zerrat(px) {
        if rat_lt(y, &rat_zero(), precision)? {
            return Err(CalcError::Domain);
        } else if zerrat(y) {
            duprat(px, &rat_one());
        }
        return Ok(());
    }

    let mut pxint = rat_zero();
    duprat(&mut pxint, px);
    _subrat(&mut pxint, &rat_one(), precision)?;

    let near_one = rat_gt(&pxint, &rat_negsmallest(), precision)?
        && rat_lt(&pxint, &rat_smallest(), precision)?
        && (sign == 1);

    if near_one {
        duprat(px, &rat_one());
        return Ok(());
    }

    // Check if y is integral
    let mut podd = rat_zero();
    duprat(&mut podd, y);
    fracrat(&mut podd, radix, precision)?;

    let effective_sign;
    if rat_gt(&podd, &rat_negsmallest(), precision)?
        && rat_lt(&podd, &rat_smallest(), precision)?
    {
        // Integer exponent: use ratpowi32
        let mut iy = rat_zero();
        duprat(&mut iy, y);
        _subrat(&mut iy, &podd, precision)?;
        let inty = rattoi32(&iy, radix, precision);

        // Domain check via log
        let mut plnx = rat_zero();
        duprat(&mut plnx, px);
        _lograt(&mut plnx, precision)?;
        mulrat(&mut plnx, &iy, precision)?;
        if rat_gt(&plnx, &rat_max_exp(), precision)?
            || rat_lt(&plnx, &rat_min_exp(), precision)?
        {
            return Err(CalcError::Domain);
        }

        ratpowi32(px, inty, precision)?;
        effective_sign = if (inty & 1) == 0 { 1 } else { sign };
    } else {
        // Fractional exponent
        effective_sign = if sign == -1 {
            // Check denominator parity for complex result
            let mut pnum = rat_zero();
            duprat(&mut pnum, &rat_zero());
            let mut pden = rat_zero();
            duprat(&mut pden, &rat_zero());
            dupnum(&mut pnum.borrow_mut().pp, &y.borrow().pp);
            pnum.borrow().pp.borrow_mut().sign = 1;
            dupnum(&mut pden.borrow_mut().pp, &y.borrow().pq);
            pden.borrow().pp.borrow_mut().sign = 1;

            let mut s = sign;
            while is_even(&pnum, radix, precision)? && is_even(&pden, radix, precision)? {
                divrat(&mut pnum, &rat_two(), precision)?;
                divrat(&mut pden, &rat_two(), precision)?;
            }
            if is_even(&pden, radix, precision)? {
                return Err(CalcError::Domain);
            }
            if is_even(&pnum, radix, precision)? {
                s = 1;
            }
            s
        } else {
            1
        };

        _lograt(px, precision)?;
        mulrat(px, y, precision)?;
        exprat(px, radix, precision)?;
    }

    px.borrow().pp.borrow_mut().sign *= effective_sign;
    Ok(())
}

// ---------------------------------------------------------------------------
//  powratNumeratorDenominator: precise power via num/denom decomposition
// ---------------------------------------------------------------------------

pub fn powrat_numerator_denominator(
    px: &mut PRat,
    y: &PRat,
    radix: u32,
    precision: i32,
) -> CalcResult<()> {
    let mut y_num = rat_zero();
    duprat(&mut y_num, &rat_zero());
    let mut y_den = rat_zero();
    duprat(&mut y_den, &rat_zero());
    dupnum(&mut y_num.borrow_mut().pp, &y.borrow().pp);
    dupnum(&mut y_den.borrow_mut().pp, &y.borrow().pq);

    let mut px_pow = rat_zero();
    duprat(&mut px_pow, px);

    if !rat_equ(&y_num, &rat_one(), precision)? {
        powratcomp(&mut px_pow, &y_num, radix, precision)?;
    }

    if !rat_equ(&y_den, &rat_one(), precision)? {
        let mut one_over_y_den = rat_one();
        divrat(&mut one_over_y_den, &y_den, precision)?;

        let mut original_result = rat_zero();
        duprat(&mut original_result, &px_pow);
        powratcomp(&mut original_result, &one_over_y_den, radix, precision)?;

        let mut rounded = rat_zero();
        duprat(&mut rounded, &original_result);
        if rounded.borrow().pp.borrow().sign == -1 {
            _subrat(&mut rounded, &rat_half(), precision)?;
        } else {
            _addrat(&mut rounded, &rat_half(), precision)?;
        }
        intrat(&mut rounded, radix, precision)?;

        let mut rounded_power = rat_zero();
        duprat(&mut rounded_power, &rounded);
        powratcomp(&mut rounded_power, &y_den, radix, precision)?;

        if rat_equ(&rounded_power, &px_pow, precision)? {
            duprat(px, &rounded);
        } else {
            duprat(px, &original_result);
        }
    } else {
        duprat(px, &px_pow);
    }

    Ok(())
}

// ---------------------------------------------------------------------------
//  powrat: top-level power dispatcher
// ---------------------------------------------------------------------------

pub fn powrat(px: &mut PRat, y: &PRat, radix: u32, precision: i32) -> CalcResult<()> {
    if zerrat(px) || zerrat(y) {
        return powratcomp(px, y, radix, precision);
    }
    if rat_equ(y, &rat_one(), precision)? {
        return Ok(());
    }

    match powrat_numerator_denominator(px, y, radix, precision) {
        Ok(()) => {}
        Err(_) => {
            powratcomp(px, y, radix, precision)?;
        }
    }
    Ok(())
}

// Need these in support::thread_local context:
use crate::ratpack::support::{rat_equ, rat_le, LN_TEN, LN_TWO};
