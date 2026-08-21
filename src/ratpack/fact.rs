// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `fact.cpp`.
//! Contains factorial and the underlying _gamma function for rationals.

use crate::ratpack::constants::{rat_max_fact, rat_min_fact, rat_one, rat_two, rat_zero};
use crate::ratpack::conv::{i32tonum, i32torat};
use crate::ratpack::errors::{CalcError, CalcResult};
use crate::ratpack::exp::{_lograt, exprat, powratcomp};
use crate::ratpack::num::{mulnumx, ratpowi32};
use crate::ratpack::rat::{_addrat, _subrat, divrat, fracrat, mulrat};
use crate::ratpack::support::{
    absrat, dupnum, duprat, logratradix, rat_gt, rat_lt, rat_neq, zerrat,
};
use crate::ratpack::types::{PRat, BASEX};
use std::rc::Rc;

// ---------------------------------------------------------------------------
//  _gamma: Lanczos-style series for the Gamma function.
//  Called with n in range (0, 1.5] after factoring out integers.
// ---------------------------------------------------------------------------

pub fn _gamma(pn: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let mut ratprec = i32torat(precision);

    // a = ln(radix) * precision + 2
    let mut a = i32torat(radix as i32);
    _lograt(&mut a, precision)?;
    mulrat(&mut a, &ratprec, precision)?;
    _addrat(&mut a, &rat_two(), precision)?;

    // a += n * ln(a) + 1
    let mut tmp = rat_zero();
    duprat(&mut tmp, &a);
    _lograt(&mut tmp, precision)?;
    mulrat(&mut tmp, pn, precision)?;
    _addrat(&mut a, &tmp, precision)?;
    _addrat(&mut a, &rat_one(), precision)?;

    // Bump precision: precision += round(ln(exp(a) * a^(n+1.5)) - ln(radix))
    {
        let mut tmp2 = rat_zero();
        duprat(&mut tmp2, pn);
        let mut one_pt_five = i32torat(3);
        divrat(&mut one_pt_five, &rat_two(), precision)?;
        _addrat(&mut tmp2, &one_pt_five, precision)?;

        let mut term = rat_zero();
        duprat(&mut term, &a);
        powratcomp(&mut term, &tmp2, radix, precision)?;

        let mut tmp3 = rat_zero();
        duprat(&mut tmp3, &a);
        exprat(&mut tmp3, radix, precision)?;
        mulrat(&mut term, &tmp3, precision)?;
        _lograt(&mut term, precision)?;

        let mut tmp4 = i32torat(radix as i32);
        _lograt(&mut tmp4, precision)?;
        _subrat(&mut term, &tmp4, precision)?;

        let extra = crate::ratpack::conv::rattoi32(&term, radix, precision);
        let precision = precision + extra;
        // (local shadow of precision now has bumped value for the rest of _gamma)
        return _gamma_inner(pn, radix, precision, a);
    }
}

fn _gamma_inner(pn: &mut PRat, radix: u32, precision: i32, a: PRat) -> CalcResult<()> {
    use crate::ratpack::conv::{i32tonum, rattoi32};
    use crate::ratpack::support::rat_gt;

    let ratprec = i32torat(precision);
    let mut factorial = rat_one();
    let mut count = i32tonum(0, BASEX);

    let mut mpy = rat_zero();
    duprat(&mut mpy, &a);
    powratcomp(&mut mpy, pn, radix, precision)?;

    // a2 = a^2
    let mut a2 = rat_zero();
    duprat(&mut a2, &a);
    mulrat(&mut a2, &a, precision)?;

    // sum = 1/n - a/(n+1)
    let mut sum = rat_one();
    divrat(&mut sum, pn, precision)?;
    let mut tmp = rat_zero();
    duprat(&mut tmp, pn);
    _addrat(&mut tmp, &rat_one(), precision)?;
    let mut term_init = rat_zero();
    duprat(&mut term_init, &a);
    divrat(&mut term_init, &tmp, precision)?;
    _subrat(&mut sum, &term_init, precision)?;

    // err = radix^(-precision) / radix
    let mut err = i32torat(radix as i32);
    let mut neg_prec = i32torat(-precision);
    powratcomp(&mut err, &neg_prec, radix, precision)?;
    divrat(&mut err, &i32torat(radix as i32), precision)?;

    let mut term = rat_two(); // Something not tiny

    while !zerrat(&term) && rat_gt(&term, &err, precision)? {
        // pn += 2
        _addrat(pn, &rat_two(), precision)?;

        // factorial: multiply by next two even numbers
        let one = i32tonum(1, BASEX);
        crate::ratpack::num::addnum(&mut count, &one, BASEX);
        mulnumx(&mut factorial.borrow_mut().pp, &count);
        let one2 = i32tonum(1, BASEX);
        crate::ratpack::num::addnum(&mut count, &one2, BASEX);
        mulnumx(&mut factorial.borrow_mut().pp, &count);
        divrat(&mut factorial, &a2, precision)?;

        let mut tmp2 = rat_zero();
        duprat(&mut tmp2, pn);
        _addrat(&mut tmp2, &rat_one(), precision)?;

        let mut new_term = rat_zero();
        // new_term = count as rat
        dupnum(&mut new_term.borrow_mut().pp, &count);
        dupnum(&mut new_term.borrow_mut().pq, &crate::ratpack::constants::num_one());
        _addrat(&mut new_term, &rat_one(), precision)?;
        mulrat(&mut new_term, &tmp2, precision)?;

        let mut tmp3 = rat_zero();
        duprat(&mut tmp3, &a);
        divrat(&mut tmp3, &new_term, precision)?;

        duprat(&mut term, &rat_one());
        divrat(&mut term, pn, precision)?;
        _subrat(&mut term, &tmp3, precision)?;
        divrat(&mut term, &factorial, precision)?;
        _addrat(&mut sum, &term, precision)?;
        absrat(&term);
    }

    mulrat(&mut sum, &mpy, precision)?;
    duprat(pn, &sum);
    Ok(())
}

// ---------------------------------------------------------------------------
//  factrat: factorial dispatcher
// ---------------------------------------------------------------------------

pub fn factrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    if rat_gt(px, &rat_max_fact(), precision)? || rat_lt(px, &rat_min_fact(), precision)? {
        return Err(CalcError::Overflow);
    }

    let mut fact = rat_one();
    let mut neg_rat_one = rat_one();
    neg_rat_one.borrow().pp.borrow_mut().sign *= -1;

    let mut frac = rat_zero();
    duprat(&mut frac, px);
    fracrat(&mut frac, radix, precision)?;

    // Negative integers are undefined
    if (zerrat(&frac) || (logratradix(&frac) <= -precision))
        && (crate::ratpack::support::rat_sign(px) == -1)
    {
        return Err(CalcError::Domain);
    }

    // Bring x down to (0, 1]
    while rat_gt(px, &rat_zero(), precision)?
        && (logratradix(px) > -precision)
    {
        mulrat(&mut fact, px, precision)?;
        _subrat(px, &rat_one(), precision)?;
    }

    // Round if very close to an integer
    if logratradix(px) <= -precision {
        duprat(px, &rat_zero());
        crate::ratpack::support::intrat(&mut fact, radix, precision)?;
    }

    // Bring x up to (-1, 0]
    while rat_lt(px, &neg_rat_one, precision)? {
        _addrat(px, &rat_one(), precision)?;
        divrat(&mut fact, px, precision)?;
    }

    if rat_neq(px, &rat_zero(), precision)? {
        _addrat(px, &rat_one(), precision)?;
        _gamma(px, radix, precision)?;
        mulrat(px, &fact, precision)?;
    } else {
        duprat(px, &fact);
    }

    Ok(())
}
