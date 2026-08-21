use crate::ratpack::types::{Number, Rat, PNumber, PRat};
use std::rc::Rc;
use std::cell::RefCell;

fn new_pnumber(sign: i32, exp: i32, mant: Vec<u32>) -> PNumber {
    Rc::new(RefCell::new(Number { sign, exp, mant }))
}

fn new_prat(pp: PNumber, pq: PNumber) -> PRat {
    Rc::new(RefCell::new(Rat { pp, pq }))
}

pub fn num_one() -> PNumber {
    new_pnumber(1, 0, vec![1])
}

pub fn num_two() -> PNumber {
    new_pnumber(1, 0, vec![2])
}

pub fn num_five() -> PNumber {
    new_pnumber(1, 0, vec![5])
}

pub fn num_six() -> PNumber {
    new_pnumber(1, 0, vec![6])
}

pub fn num_ten() -> PNumber {
    new_pnumber(1, 0, vec![10])
}

pub fn rat_zero() -> PRat {
    new_prat(new_pnumber(1, 0, vec![0]), new_pnumber(1, 0, vec![1]))
}

pub fn rat_one() -> PRat {
    new_prat(new_pnumber(1, 0, vec![1]), new_pnumber(1, 0, vec![1]))
}

pub fn rat_two() -> PRat {
    new_prat(new_pnumber(1, 0, vec![2]), new_pnumber(1, 0, vec![1]))
}

pub fn rat_half() -> PRat {
    new_prat(new_pnumber(1, 0, vec![1]), new_pnumber(1, 0, vec![2]))
}

pub fn rat_neg_one() -> PRat {
    new_prat(new_pnumber(-1, 0, vec![1]), new_pnumber(1, 0, vec![1]))
}
