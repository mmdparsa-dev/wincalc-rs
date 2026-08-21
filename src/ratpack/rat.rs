use crate::ratpack::types::*;
use crate::ratpack::constants::*;
use crate::ratpack::*;
use std::rc::Rc;
use std::cell::RefCell;

pub fn renormalize(x: &PRat) {
    let mut pp = x.borrow().pp.borrow_mut();
    let mut pq = x.borrow().pq.borrow_mut();
    if pp.exp < 0 {
        pq.exp -= pp.exp;
        pp.exp = 0;
    }
    if pq.exp < 0 {
        pp.exp -= pq.exp;
        pq.exp = 0;
    }
}

pub fn dupnum(num: &PNumber) -> PNumber {
    Rc::new(RefCell::new(num.borrow().clone()))
}

pub fn duprat(rat: &PRat) -> PRat {
    Rc::new(RefCell::new(Rat {
        pp: dupnum(&rat.borrow().pp),
        pq: dupnum(&rat.borrow().pq),
    }))
}

pub fn absrat(x: &PRat) {
    x.borrow().pp.borrow_mut().sign = 1;
    x.borrow().pq.borrow_mut().sign = 1;
}

pub fn gcdrat(pa: &mut PRat, precision: i32) {
    let a = Rc::clone(pa);
    let pgcd = gcd(&a.borrow().pp, &a.borrow().pq);

    if !zernum(&pgcd) {
        let mut pp = Rc::clone(&a.borrow().pp);
        divnumx(&mut pp, &pgcd, precision);
        a.borrow_mut().pp = pp;

        let mut pq = Rc::clone(&a.borrow().pq);
        divnumx(&mut pq, &pgcd, precision);
        a.borrow_mut().pq = pq;
    }

    *pa = a;
    renormalize(pa);
}

pub fn fracrat(pa: &mut PRat, radix: u32, precision: i32) {
    if !zernum(&pa.borrow().pp) && !equnum(&pa.borrow().pq, &num_one()) {
        let mut a = Rc::clone(pa);
        flatrat(&mut a, radix, precision);
        *pa = a;
    }

    let mut pp = Rc::clone(&pa.borrow().pp);
    remnum(&mut pp, &pa.borrow().pq, BASEX);
    pa.borrow_mut().pp = pp;

    renormalize(pa);
}

pub fn mulrat(pa: &mut PRat, b: &PRat, precision: i32) {
    if !zernum(&pa.borrow().pp) {
        let mut pp = Rc::clone(&pa.borrow().pp);
        mulnumx(&mut pp, &b.borrow().pp);
        pa.borrow_mut().pp = pp;

        let mut pq = Rc::clone(&pa.borrow().pq);
        mulnumx(&mut pq, &b.borrow().pq);
        pa.borrow_mut().pq = pq;

        trimit(pa, precision);
    } else {
        pa.borrow_mut().pq = dupnum(&num_one());
    }

    // #[cfg(feature = "mulgcd")]
    // gcdrat(pa, precision);
}

pub fn divrat(pa: &mut PRat, b: &PRat, precision: i32) {
    if !zernum(&pa.borrow().pp) {
        let mut pp = Rc::clone(&pa.borrow().pp);
        mulnumx(&mut pp, &b.borrow().pq);
        pa.borrow_mut().pp = pp;

        let mut pq = Rc::clone(&pa.borrow().pq);
        mulnumx(&mut pq, &b.borrow().pp);
        pa.borrow_mut().pq = pq;

        if zernum(&pa.borrow().pq) {
            panic!("CALC_E_DIVIDEBYZERO");
        }
        trimit(pa, precision);
    } else {
        if zerrat(b) {
            panic!("CALC_E_INDEFINITE");
        } else {
            pa.borrow_mut().pq = dupnum(&num_one());
        }
    }

    // #[cfg(feature = "divgcd")]
    // gcdrat(pa, precision);
}

pub fn subrat(pa: &mut PRat, b: &PRat, precision: i32) {
    let a = duprat(pa);
    _subrat(pa, b, precision);
    _snaprat(pa, &a, Some(b), precision);
}

pub fn _subrat(pa: &mut PRat, b: &PRat, precision: i32) {
    b.borrow().pp.borrow_mut().sign *= -1;
    _addrat(pa, b, precision);
    b.borrow().pp.borrow_mut().sign *= -1;
}

pub fn addrat(pa: &mut PRat, b: &PRat, precision: i32) {
    let a = duprat(pa);
    _addrat(pa, b, precision);
    _snaprat(pa, &a, Some(b), precision);
}

pub fn _addrat(pa: &mut PRat, b: &PRat, precision: i32) {
    let q_eq = equnum(&pa.borrow().pq, &b.borrow().pq);

    if q_eq {
        let pa_pq_sign = pa.borrow().pq.borrow().sign;
        pa.borrow().pp.borrow_mut().sign *= pa_pq_sign;
        pa.borrow().pq.borrow_mut().sign = 1;

        let b_pq_sign = b.borrow().pq.borrow().sign;
        b.borrow().pp.borrow_mut().sign *= b_pq_sign;
        b.borrow().pq.borrow_mut().sign = 1;

        let mut pp = Rc::clone(&pa.borrow().pp);
        addnum(&mut pp, &b.borrow().pp, BASEX);
        pa.borrow_mut().pp = pp;
    } else {
        let mut bot = dupnum(&pa.borrow().pq);
        mulnumx(&mut bot, &b.borrow().pq);

        let mut pp = Rc::clone(&pa.borrow().pp);
        mulnumx(&mut pp, &b.borrow().pq);
        pa.borrow_mut().pp = pp;

        let mut pq = Rc::clone(&pa.borrow().pq);
        mulnumx(&mut pq, &b.borrow().pp);
        pa.borrow_mut().pq = pq;

        let mut pp_again = Rc::clone(&pa.borrow().pp);
        addnum(&mut pp_again, &pa.borrow().pq, BASEX);
        pa.borrow_mut().pp = pp_again;

        pa.borrow_mut().pq = bot;
        trimit(pa, precision);

        let pq_sign = pa.borrow().pq.borrow().sign;
        pa.borrow().pp.borrow_mut().sign *= pq_sign;
        pa.borrow().pq.borrow_mut().sign = 1;
    }

    // #[cfg(feature = "addgcd")]
    // gcdrat(pa, precision);
}

pub fn rootrat(py: &mut PRat, n: &PRat, radix: u32, precision: i32) {
    let mut oneovern = duprat(&rat_one());
    divrat(&mut oneovern, n, precision);
    powrat(py, &oneovern, radix, precision);
}

pub fn zerrat(a: &PRat) -> bool {
    zernum(&a.borrow().pp)
}

pub fn _snaprat(pr: &mut PRat, a: &PRat, b: Option<&PRat>, precision: i32) {
    let mut threshold = if let Some(b_val) = b {
        let abs_a = duprat(a);
        let abs_b = duprat(b_val);
        absrat(&abs_a);
        absrat(&abs_b);

        if rat_lt(&abs_a, &abs_b, precision) {
            duprat(&abs_b)
        } else {
            duprat(&abs_a)
        }
    } else {
        let t = duprat(a);
        absrat(&t);
        t
    };

    mulrat(&mut threshold, &rat_smallest(), precision);

    let abs_r = duprat(pr);
    absrat(&abs_r);

    if rat_lt(&abs_r, &threshold, precision) {
        *pr = duprat(&rat_zero());
    }
}
