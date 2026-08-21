// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `transh.cpp`.
//! Contains hyperbolic sin, cos, and tan for rationals.

use crate::ratpack::constants::{rat_min_exp, rat_one, rat_ten, rat_two, rat_zero};
use crate::ratpack::errors::{CalcError, CalcResult};
use crate::ratpack::exp::exprat;
use crate::ratpack::num::addnum;
use crate::ratpack::rat::{_addrat, _subrat, divrat};
use crate::ratpack::support::{dupnum, duprat, rat_ge, rat_gt, rat_lt, trimit, zerrat};
use crate::ratpack::types::{PRat, BASEX};
use crate::ratpack::conv::i32tonum;
use std::rc::Rc;

/// IsValidForHypFunc: checks that x > rat_min_exp / 10.
fn is_valid_for_hyp_func(px: &PRat, precision: i32) -> CalcResult<bool> {
    let mut ptmp = rat_zero();
    duprat(&mut ptmp, &rat_min_exp());
    divrat(&mut ptmp, &rat_ten(), precision)?;
    if rat_lt(px, &ptmp, precision)? {
        return Ok(false);
    }
    Ok(true)
}

// ---------------------------------------------------------------------------
//  _sinhrat: Taylor series for sinh(x)
// ---------------------------------------------------------------------------

pub fn _sinhrat(px: &mut PRat, precision: i32) -> CalcResult<()> {
    if !is_valid_for_hyp_func(px, precision)? {
        return Err(CalcError::Domain);
    }

    let mut pret = rat_zero();
    duprat(&mut pret, px);
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, &pret);

    // xx = x^2 (positive, unlike sin which negates)
    let mut xx = rat_zero();
    duprat(&mut xx, px);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pp, &px.borrow().pp);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pq, &px.borrow().pq);

    let mut n2 = i32tonum(1, BASEX);

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
    Ok(())
}

pub fn sinhrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    if rat_ge(px, &rat_one(), precision)? {
        let mut tmpx = rat_zero();
        duprat(&mut tmpx, px);
        exprat(px, radix, precision)?;
        tmpx.borrow().pp.borrow_mut().sign *= -1;
        exprat(&mut tmpx, radix, precision)?;
        _subrat(px, &tmpx, precision)?;
        divrat(px, &rat_two(), precision)?;
    } else {
        _sinhrat(px, precision)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
//  _coshrat: Taylor series for cosh(x)
// ---------------------------------------------------------------------------

pub fn _coshrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    if !is_valid_for_hyp_func(px, precision)? {
        return Err(CalcError::Domain);
    }

    let mut pret = crate::ratpack::conv::i32torat(1);
    let mut thisterm = rat_zero();
    duprat(&mut thisterm, &pret);

    let mut xx = rat_zero();
    duprat(&mut xx, px);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pp, &px.borrow().pp);
    crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pq, &px.borrow().pq);

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
    Ok(())
}

pub fn coshrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    // cosh is symmetric: take absolute value first
    px.borrow().pp.borrow_mut().sign = 1;
    px.borrow().pq.borrow_mut().sign = 1;

    if rat_ge(px, &rat_one(), precision)? {
        let mut tmpx = rat_zero();
        duprat(&mut tmpx, px);
        exprat(px, radix, precision)?;
        tmpx.borrow().pp.borrow_mut().sign *= -1;
        exprat(&mut tmpx, radix, precision)?;
        _addrat(px, &tmpx, precision)?;
        divrat(px, &rat_two(), precision)?;
    } else {
        _coshrat(px, radix, precision)?;
    }

    // Snap to >= 1 (cosh >= 1 always)
    if crate::ratpack::support::rat_lt(px, &rat_one(), precision)? {
        duprat(px, &rat_one());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
//  tanhrat: tanh = sinh / cosh
// ---------------------------------------------------------------------------

pub fn tanhrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let mut ptmp = rat_zero();
    duprat(&mut ptmp, px);
    sinhrat(px, radix, precision)?;
    coshrat(&mut ptmp, radix, precision)?;

    // Cross-multiply instead of full divrat (mirrors C++)
    crate::ratpack::num::mulnumx(&mut px.borrow_mut().pp, &ptmp.borrow().pq);
    crate::ratpack::num::mulnumx(&mut px.borrow_mut().pq, &ptmp.borrow().pp);

    Ok(())
}
