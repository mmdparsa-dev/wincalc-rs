// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `trans.cpp`.
//! Contains sin, cos, and tan for rationals.

use crate::ratpack::constants::{rat_negsmallest, rat_one, rat_smallest, rat_zero};
use crate::ratpack::errors::CalcResult;
use crate::ratpack::exp::{_exprat, exprat};
use crate::ratpack::rat::{_addrat, _subrat, divrat, mulrat};
use crate::ratpack::support::{
    dupnum, duprat, inbetween, rat_ge, rat_gt, rat_le, rat_lt, scale, scale2pi, trimit,
};
use crate::ratpack::types::{AngleType, PRat, BASEX};
use crate::ratpack::conv::i32tonum;
use crate::ratpack::num::addnum;
use std::rc::Rc;

// Globals accessed as thread_local in support.rs
macro_rules! get_const {
    ($name:ident) => {
        crate::ratpack::support::$name.with(|v| Rc::clone(&v.borrow()))
    };
}

pub fn scalerat(
    pa: &mut PRat,
    angletype: AngleType,
    radix: u32,
    precision: i32,
) -> CalcResult<()> {
    match angletype {
        AngleType::Radians => scale2pi(pa, radix, precision),
        AngleType::Degrees => {
            let r360 = crate::ratpack::support::RAT_360
                .with(|v| Rc::clone(&v.borrow()));
            scale(pa, &r360, radix, precision)
        }
        AngleType::Gradians => {
            let r400 = crate::ratpack::support::RAT_400
                .with(|v| Rc::clone(&v.borrow()));
            scale(pa, &r400, radix, precision)
        }
    }
}

// ---------------------------------------------------------------------------
//  _sinrat: Taylor series for sin(x) starting from x already in [0, 2π)
// ---------------------------------------------------------------------------

pub fn _sinrat(px: &mut PRat, precision: i32) -> CalcResult<()> {
    let mut pret = rat_zero();
    duprat(&mut pret, px);
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, px);

    // xx = -x^2 (for the denominator increments)
    let mut xx = rat_zero();
    duprat(&mut xx, px);
    dupnum(
        &mut xx.borrow_mut().pp,
        &{
            let old = Rc::clone(&xx.borrow().pp);
            let mut n = old.borrow().clone();
            n.sign *= -1;
            Rc::new(std::cell::RefCell::new(n))
        },
    );
    // Actually xx = x^2 negated for sin series: multiply x by itself then negate
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pp, &px.borrow().pp);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pq, &px.borrow().pq);
    xx.borrow().pp.borrow_mut().sign *= -1;

    let mut n2 = i32tonum(1, BASEX);

    loop {
        // NEXTTERM: thisterm *= xx / (n2*(n2+1))
        let one = i32tonum(1, BASEX);
        addnum(&mut n2, &one, BASEX); // n2++
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &xx.borrow().pp);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &xx.borrow().pq);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &n2);
        let one2 = i32tonum(1, BASEX);
        addnum(&mut n2, &one2, BASEX); // n2++
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &n2);

        trimit(&mut thisterm, precision)?;
        _addrat(&mut pret, &thisterm, precision)?;

        if crate::ratpack::support::lograt2(&thisterm) < -precision {
            break;
        }
    }

    duprat(px, &pret);

    // Snap to [-1, 1]
    inbetween(px, &rat_one(), precision)?;
    // Snap near-zero to exactly zero
    if rat_le(px, &rat_smallest(), precision)? && rat_ge(px, &rat_negsmallest(), precision)? {
        duprat(px, &rat_zero());
    }
    Ok(())
}

pub fn sinrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    scale2pi(px, radix, precision)?;
    _sinrat(px, precision)
}

pub fn sinanglerat(pa: &mut PRat, angletype: AngleType, radix: u32, precision: i32) -> CalcResult<()> {
    scalerat(pa, angletype, radix, precision)?;
    match angletype {
        AngleType::Degrees => {
            let r180 = crate::ratpack::support::RAT_180.with(|v| Rc::clone(&v.borrow()));
            let r360 = crate::ratpack::support::RAT_360.with(|v| Rc::clone(&v.borrow()));
            let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
            if rat_gt(pa, &r180, precision)? {
                _subrat(pa, &r360, precision)?;
            }
            divrat(pa, &r180, precision)?;
            mulrat(pa, &pi, precision)?;
        }
        AngleType::Gradians => {
            let r200 = crate::ratpack::support::RAT_200.with(|v| Rc::clone(&v.borrow()));
            let r400 = crate::ratpack::support::RAT_400.with(|v| Rc::clone(&v.borrow()));
            let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
            if rat_gt(pa, &r200, precision)? {
                _subrat(pa, &r400, precision)?;
            }
            divrat(pa, &r200, precision)?;
            mulrat(pa, &pi, precision)?;
        }
        AngleType::Radians => {}
    }
    _sinrat(pa, precision)
}

// ---------------------------------------------------------------------------
//  _cosrat: Taylor series for cos(x)
// ---------------------------------------------------------------------------

pub fn _cosrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    use crate::ratpack::conv::i32tonum;

    let mut pret = crate::ratpack::conv::i32torat(1);
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, &pret);

    let mut xx = rat_zero();
    duprat(&mut xx, px);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pp, &px.borrow().pp);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pq, &px.borrow().pq);
    xx.borrow().pp.borrow_mut().sign *= -1;

    let mut n2 = i32tonum(0, BASEX);

    loop {
        let one = i32tonum(1, BASEX);
        addnum(&mut n2, &one, BASEX);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pp, &xx.borrow().pp);
        crate::ratpack::num::mulnumx(&mut thisterm.borrow_mut().pq, &xx.borrow().pq);
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

    inbetween(px, &rat_one(), precision)?;
    if rat_le(px, &rat_smallest(), precision)? && rat_ge(px, &rat_negsmallest(), precision)? {
        duprat(px, &rat_zero());
    }
    Ok(())
}

pub fn cosrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    scale2pi(px, radix, precision)?;
    _cosrat(px, radix, precision)
}

pub fn cosanglerat(pa: &mut PRat, angletype: AngleType, radix: u32, precision: i32) -> CalcResult<()> {
    scalerat(pa, angletype, radix, precision)?;
    match angletype {
        AngleType::Degrees => {
            let r180 = crate::ratpack::support::RAT_180.with(|v| Rc::clone(&v.borrow()));
            let r360 = crate::ratpack::support::RAT_360.with(|v| Rc::clone(&v.borrow()));
            let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
            if rat_gt(pa, &r180, precision)? {
                let mut ptmp = rat_zero();
                duprat(&mut ptmp, &r360);
                _subrat(&mut ptmp, pa, precision)?;
                duprat(pa, &ptmp);
            }
            divrat(pa, &r180, precision)?;
            mulrat(pa, &pi, precision)?;
        }
        AngleType::Gradians => {
            let r200 = crate::ratpack::support::RAT_200.with(|v| Rc::clone(&v.borrow()));
            let r400 = crate::ratpack::support::RAT_400.with(|v| Rc::clone(&v.borrow()));
            let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
            if rat_gt(pa, &r200, precision)? {
                let mut ptmp = rat_zero();
                duprat(&mut ptmp, &r400);
                _subrat(&mut ptmp, pa, precision)?;
                duprat(pa, &ptmp);
            }
            divrat(pa, &r200, precision)?;
            mulrat(pa, &pi, precision)?;
        }
        AngleType::Radians => {}
    }
    _cosrat(pa, radix, precision)
}

// ---------------------------------------------------------------------------
//  _tanrat / tanrat: tan = sin/cos
// ---------------------------------------------------------------------------

pub fn _tanrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let mut ptmp = rat_zero();
    duprat(&mut ptmp, px);
    _sinrat(px, precision)?;
    _cosrat(&mut ptmp, radix, precision)?;
    if crate::ratpack::support::zerrat(&ptmp) {
        return Err(crate::ratpack::errors::CalcError::Domain);
    }
    crate::ratpack::rat::divrat(px, &ptmp, precision)?;
    Ok(())
}

pub fn tanrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    scale2pi(px, radix, precision)?;
    _tanrat(px, radix, precision)
}

pub fn tananglerat(pa: &mut PRat, angletype: AngleType, radix: u32, precision: i32) -> CalcResult<()> {
    scalerat(pa, angletype, radix, precision)?;
    match angletype {
        AngleType::Degrees => {
            let r180 = crate::ratpack::support::RAT_180.with(|v| Rc::clone(&v.borrow()));
            let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
            if rat_gt(pa, &r180, precision)? {
                _subrat(pa, &r180, precision)?;
            }
            divrat(pa, &r180, precision)?;
            mulrat(pa, &pi, precision)?;
        }
        AngleType::Gradians => {
            let r200 = crate::ratpack::support::RAT_200.with(|v| Rc::clone(&v.borrow()));
            let pi = crate::ratpack::support::PI.with(|v| Rc::clone(&v.borrow()));
            if rat_gt(pa, &r200, precision)? {
                _subrat(pa, &r200, precision)?;
            }
            divrat(pa, &r200, precision)?;
            mulrat(pa, &pi, precision)?;
        }
        AngleType::Radians => {}
    }
    _tanrat(pa, radix, precision)
}
