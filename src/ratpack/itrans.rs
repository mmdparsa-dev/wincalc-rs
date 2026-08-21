// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `itrans.cpp`.
//! Contains inverse sin, cos, and tan for rationals.

use crate::ratpack::constants::{pi_over_two, pt_eight_five, rat_negsmallest, rat_one, rat_smallest, rat_two, rat_zero};
use crate::ratpack::errors::{CalcError, CalcResult};
use crate::ratpack::exp::rootrat;
use crate::ratpack::num::addnum;
use crate::ratpack::rat::{_addrat, _subrat, divrat, mulrat};
use crate::ratpack::support::{dupnum, duprat, rat_equ, rat_ge, rat_gt, rat_le, rat_lt, trimit};
use crate::ratpack::types::{AngleType, PRat, BASEX};
use crate::ratpack::conv::i32tonum;
use std::rc::Rc;

pub fn ascalerat(pa: &mut PRat, angletype: AngleType, precision: i32) -> CalcResult<()> {
    let two_pi = crate::ratpack::support::TWO_PI.with(|v| Rc::clone(&v.borrow()));
    match angletype {
        AngleType::Radians => {}
        AngleType::Degrees => {
            divrat(pa, &two_pi, precision)?;
            let r360 = crate::ratpack::support::RAT_360.with(|v| Rc::clone(&v.borrow()));
            mulrat(pa, &r360, precision)?;
        }
        AngleType::Gradians => {
            divrat(pa, &two_pi, precision)?;
            let r400 = crate::ratpack::support::RAT_400.with(|v| Rc::clone(&v.borrow()));
            mulrat(pa, &r400, precision)?;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
//  _asinrat: Taylor series for asin(x) for |x| <= 0.85
// ---------------------------------------------------------------------------

pub fn _asinrat(px: &mut PRat, precision: i32) -> CalcResult<()> {
    let mut pret = rat_zero();
    duprat(&mut pret, px);
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, px);
    let mut n2 = i32tonum(1, BASEX);

    // xx = x^2 (positive, ratio multiplied per step)
    let mut xx = rat_zero();
    duprat(&mut xx, px);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pp, &px.borrow().pp);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pq, &px.borrow().pq);

    loop {
        // NEXTTERM: thisterm *= n2^2 * x^2 / ((n2+1)*(n2+2))
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &n2);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &n2);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &xx.borrow().pp);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &xx.borrow().pq);
        let one = i32tonum(1, BASEX);
        addnum(&mut n2, &one, BASEX);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &n2);
        let one2 = i32tonum(1, BASEX);
        addnum(&mut n2, &one2, BASEX);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &n2);

        trimit(&mut thisterm, precision)?;
        _addrat(&mut pret, &thisterm, precision)?;

        if crate::ratpack::support::lograt2(&thisterm) < -precision {
            break;
        }
    }

    duprat(px, &pret);
    Ok(())
}

pub fn asinrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let sgn = crate::ratpack::support::rat_sign(px);
    px.borrow().pp.borrow_mut().sign = 1;
    px.borrow().pq.borrow_mut().sign = 1;

    // Check if x ≈ 1 (near boundary)
    let mut phack = rat_zero();
    duprat(&mut phack, px);
    _subrat(&mut phack, &rat_one(), precision)?;
    let near_one = rat_le(&phack, &rat_smallest(), precision)?
        && rat_ge(&phack, &rat_negsmallest(), precision)?;

    if near_one {
        duprat(px, &pi_over_two());
    } else if rat_gt(px, &pt_eight_five(), precision)? {
        // Near ±1 but not exactly: use asin(sqrt(1-x^2)) alternative form
        if rat_gt(px, &rat_one(), precision)? {
            _subrat(px, &rat_one(), precision)?;
            if rat_gt(px, &rat_smallest(), precision)? {
                return Err(CalcError::Domain);
            } else {
                duprat(px, &rat_one());
            }
        }
        let mut pret = rat_zero();
        duprat(&mut pret, px);
        mulrat(px, &pret, precision)?;
        px.borrow().pp.borrow_mut().sign *= -1;
        _addrat(px, &rat_one(), precision)?;
        rootrat(px, &rat_two(), radix, precision)?;
        _asinrat(px, precision)?;
        px.borrow().pp.borrow_mut().sign *= -1;
        let pio2 = pi_over_two();
        _addrat(px, &pio2, precision)?;
    } else {
        _asinrat(px, precision)?;
    }

    px.borrow().pp.borrow_mut().sign = sgn;
    px.borrow().pq.borrow_mut().sign = 1;
    Ok(())
}

pub fn asinanglerat(pa: &mut PRat, angletype: AngleType, radix: u32, precision: i32) -> CalcResult<()> {
    asinrat(pa, radix, precision)?;
    ascalerat(pa, angletype, precision)
}

// ---------------------------------------------------------------------------
//  acosrat: acos(x) = π/2 - asin(x)
// ---------------------------------------------------------------------------

pub fn acosrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let sgn = crate::ratpack::support::rat_sign(px);
    px.borrow().pp.borrow_mut().sign = 1;
    px.borrow().pq.borrow_mut().sign = 1;

    if rat_equ(px, &rat_one(), precision)? {
        if sgn == -1 {
            let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
            duprat(px, &pi);
        } else {
            duprat(px, &rat_zero());
        }
    } else {
        px.borrow().pp.borrow_mut().sign = sgn;
        asinrat(px, radix, precision)?;
        px.borrow().pp.borrow_mut().sign *= -1;
        _addrat(px, &pi_over_two(), precision)?;
    }
    Ok(())
}

pub fn acosanglerat(pa: &mut PRat, angletype: AngleType, radix: u32, precision: i32) -> CalcResult<()> {
    acosrat(pa, radix, precision)?;
    ascalerat(pa, angletype, precision)
}

// ---------------------------------------------------------------------------
//  _atanrat: Taylor series for atan(x) for |x| <= 0.85
// ---------------------------------------------------------------------------

pub fn _atanrat(px: &mut PRat, precision: i32) -> CalcResult<()> {
    let mut pret = rat_zero();
    duprat(&mut pret, px);
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, px);
    let mut n2 = i32tonum(1, BASEX);

    // xx = -x^2
    let mut xx = rat_zero();
    duprat(&mut xx, px);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pp, &px.borrow().pp);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pq, &px.borrow().pq);
    xx.borrow().pp.borrow_mut().sign *= -1;

    loop {
        // NEXTTERM: thisterm *= n2 * xx / (n2+2)
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &n2);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &xx.borrow().pp);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &xx.borrow().pq);
        let one = i32tonum(1, BASEX);
        addnum(&mut n2, &one, BASEX);
        let one2 = i32tonum(1, BASEX);
        addnum(&mut n2, &one2, BASEX);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &n2);

        trimit(&mut thisterm, precision)?;
        _addrat(&mut pret, &thisterm, precision)?;

        if crate::ratpack::support::lograt2(&thisterm) < -precision {
            break;
        }
    }

    duprat(px, &pret);
    Ok(())
}

pub fn atanrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let sgn = crate::ratpack::support::rat_sign(px);
    px.borrow().pp.borrow_mut().sign = 1;
    px.borrow().pq.borrow_mut().sign = 1;

    if rat_gt(px, &pt_eight_five(), precision)? {
        if rat_gt(px, &rat_two(), precision)? {
            // atan(x) = π/2 - atan(1/x) for |x| > 2
            px.borrow().pp.borrow_mut().sign = sgn;
            px.borrow().pq.borrow_mut().sign = 1;
            let mut tmpx = rat_one();
            divrat(&mut tmpx, px, precision)?;
            _atanrat(&mut tmpx, precision)?;
            tmpx.borrow().pp.borrow_mut().sign = sgn;
            tmpx.borrow().pq.borrow_mut().sign = 1;
            duprat(px, &pi_over_two());
            _subrat(px, &tmpx, precision)?;
        } else {
            // atan(x) = asin(x/sqrt(1+x^2)) for 0.85 < |x| <= 2
            px.borrow().pp.borrow_mut().sign = sgn;
            let mut tmpx = rat_zero();
            duprat(&mut tmpx, px);
            mulrat(&mut tmpx, px, precision)?;
            _addrat(&mut tmpx, &rat_one(), precision)?;
            rootrat(&mut tmpx, &rat_two(), radix, precision)?;
            divrat(px, &tmpx, precision)?;
            asinrat(px, radix, precision)?;
            px.borrow().pp.borrow_mut().sign = sgn;
            px.borrow().pq.borrow_mut().sign = 1;
        }
    } else {
        px.borrow().pp.borrow_mut().sign = sgn;
        px.borrow().pq.borrow_mut().sign = 1;
        _atanrat(px, precision)?;
    }

    let pio2 = pi_over_two();
    if rat_gt(px, &pio2, precision)? {
        let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
        _subrat(px, &pi, precision)?;
    }
    Ok(())
}

pub fn atananglerat(pa: &mut PRat, angletype: AngleType, radix: u32, precision: i32) -> CalcResult<()> {
    atanrat(pa, radix, precision)?;
    ascalerat(pa, angletype, precision)
}
