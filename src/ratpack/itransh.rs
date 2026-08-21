// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

//! Rust translation of `itransh.cpp`.
//! Contains inverse hyperbolic sin, cos, and tan for rationals.

use crate::ratpack::constants::{pt_eight_five, rat_one, rat_two, rat_zero};
use crate::ratpack::errors::{CalcError, CalcResult};
use crate::ratpack::exp::{_lograt, rootrat};
use crate::ratpack::num::addnum;
use crate::ratpack::rat::{_addrat, _subrat, divrat, mulrat};
use crate::ratpack::support::{duprat, rat_gt, rat_lt, trimit};
use crate::ratpack::types::{PRat, BASEX};
use crate::ratpack::conv::i32tonum;
use std::rc::Rc;

// ---------------------------------------------------------------------------
//  asinhrat: asinh(x) = log(x + sqrt(x^2+1)) for |x| >= 0.85
//            or Taylor series for |x| < 0.85
// ---------------------------------------------------------------------------

pub fn asinhrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    let mut neg_pt_eight_five = rat_zero();
    duprat(&mut neg_pt_eight_five, &pt_eight_five());
    neg_pt_eight_five.borrow().pp.borrow_mut().sign *= -1;

    if rat_gt(px, &pt_eight_five(), precision)?
        || rat_lt(px, &neg_pt_eight_five, precision)?
    {
        // asinh(x) = log(x + sqrt(x^2+1))
        let mut ptmp = rat_zero();
        duprat(&mut ptmp, px);
        mulrat(&mut ptmp, px, precision)?;
        _addrat(&mut ptmp, &rat_one(), precision)?;
        rootrat(&mut ptmp, &rat_two(), radix, precision)?;
        _addrat(px, &ptmp, precision)?;
        _lograt(px, precision)?;
    } else {
        // Taylor series: thisterm *= -n2^2 * x^2 / ((n2+1)*(n2+2))
        let mut pret = rat_zero();
        duprat(&mut pret, px);
        let mut thisterm = rat_zero();
        duprat(&mut thisterm, px);

        // xx = -x^2
        let mut xx = rat_zero();
        duprat(&mut xx, px);
        crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pp, &px.borrow().pp);
        crate::ratpack::num::mulnumx(&mut xx.borrow_mut().pq, &px.borrow().pq);
        xx.borrow().pp.borrow_mut().sign *= -1;

        let mut n2 = i32tonum(1, BASEX);

        loop {
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
    }
    Ok(())
}

// ---------------------------------------------------------------------------
//  acoshrat: acosh(x) = ln(x + sqrt(x^2 - 1)) for x >= 1
// ---------------------------------------------------------------------------

pub fn acoshrat(px: &mut PRat, radix: u32, precision: i32) -> CalcResult<()> {
    if rat_lt(px, &rat_one(), precision)? {
        return Err(CalcError::Domain);
    }

    let mut ptmp = rat_zero();
    duprat(&mut ptmp, px);
    mulrat(&mut ptmp, px, precision)?;
    _subrat(&mut ptmp, &rat_one(), precision)?;
    rootrat(&mut ptmp, &rat_two(), radix, precision)?;
    _addrat(px, &ptmp, precision)?;
    _lograt(px, precision)?;
    Ok(())
}

// ---------------------------------------------------------------------------
//  atanhrat: atanh(x) = (1/2) * ln((1+x)/(1-x))
// ---------------------------------------------------------------------------

pub fn atanhrat(px: &mut PRat, precision: i32) -> CalcResult<()> {
    let mut ptmp = rat_zero();
    duprat(&mut ptmp, px);
    _subrat(&mut ptmp, &rat_one(), precision)?;  // ptmp = x - 1
    _addrat(px, &rat_one(), precision)?;         // px = x + 1
    divrat(px, &ptmp, precision)?;               // px = (x+1)/(x-1)
    px.borrow().pp.borrow_mut().sign *= -1;      // px = (1+x)/(1-x)  [negate]
    _lograt(px, precision)?;
    divrat(px, &rat_two(), precision)?;
    Ok(())
}
