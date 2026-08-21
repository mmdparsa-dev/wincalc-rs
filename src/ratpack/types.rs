use std::rc::Rc;
use std::cell::RefCell;

pub const BASEXPWR: u32 = 31;
pub const BASEX: u32 = 0x80000000;
pub const MAX_LONG_SIZE: u32 = 33;

pub type MantType = u32;
pub type TwoMantType = u64;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NumberFormat {
    Float,
    Scientific,
    Engineering,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AngleType {
    Degrees,
    Radians,
    Gradians,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Number {
    pub sign: i32,
    pub exp: i32,
    pub mant: Vec<MantType>,
}

pub type PNumber = Rc<RefCell<Number>>;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Rat {
    pub pp: PNumber,
    pub pq: PNumber,
}

pub type PRat = Rc<RefCell<Rat>>;

impl Number {
    pub fn new() -> PNumber {
        Rc::new(RefCell::new(Number {
            sign: 1,
            exp: 0,
            mant: Vec::new(),
        }))
    }
}

impl Rat {
    pub fn new() -> PRat {
        Rc::new(RefCell::new(Rat {
            pp: Number::new(),
            pq: Number::new(),
        }))
    }
}
