use crate::ratpack::*;
use crate::ratpack::types::*;

pub fn mulnumx(pa: &mut PNumber, b: PNumber) {
    let b_is_not_one = {
        let b_ref = b.borrow();
        b_ref.cdigit > 1 || b_ref.mant[0] != 1 || b_ref.exp != 0
    };

    if b_is_not_one {
        let a_is_not_one = {
            let a_ref = pa.borrow();
            a_ref.cdigit > 1 || a_ref.mant[0] != 1 || a_ref.exp != 0
        };

        if a_is_not_one {
            _mulnumx(pa, b);
        } else {
            let sign = pa.borrow().sign;
            *pa = dupnum(&b);
            pa.borrow_mut().sign *= sign;
        }
    } else {
        let b_sign = b.borrow().sign;
        pa.borrow_mut().sign *= b_sign;
    }
}

pub fn _mulnumx(pa: &mut PNumber, b: PNumber) {
    let a = std::rc::Rc::clone(pa);

    let (a_cdigit, a_exp, a_sign) = {
        let a_ref = a.borrow();
        (a_ref.cdigit, a_ref.exp, a_ref.sign)
    };
    let (b_cdigit, b_exp, b_sign) = {
        let b_ref = b.borrow();
        (b_ref.cdigit, b_ref.exp, b_ref.sign)
    };

    let ibdigit = a_cdigit + b_cdigit - 1;
    let c = createnum(ibdigit + 1);

    {
        let a_ref = a.borrow();
        let b_ref = b.borrow();
        let mut c_ref = c.borrow_mut();

        c_ref.cdigit = ibdigit;
        c_ref.sign = a_sign * b_sign;
        c_ref.exp = a_exp + b_exp;

        let mut c_offset = 0;
        let mut a_idx = 0;

        for iadigit in (1..=a_cdigit).rev() {
            let da = a_ref.mant[a_idx as usize];
            a_idx += 1;
            let mut b_idx = 0;
            let mut c_idx = c_offset;
            c_offset += 1;

            for ibdigit_inner in (1..=b_cdigit).rev() {
                let mut cy: u64 = 0;
                let mut mcy: u64 = (da as u64) * (b_ref.mant[b_idx as usize] as u64);
                if mcy != 0 {
                    if ibdigit_inner == 1 && iadigit == 1 {
                        c_ref.cdigit += 1;
                    }
                }

                let mut icdigit = c_idx;
                while mcy != 0 || cy != 0 {
                    cy += (c_ref.mant[icdigit as usize] as u64) + ((mcy as u32) & (!BASEX)) as u64;
                    c_ref.mant[icdigit as usize] = (cy as u32) & (!BASEX);
                    icdigit += 1;
                    mcy >>= BASEXPWR;
                    cy >>= BASEXPWR;
                }
                b_idx += 1;
                c_idx += 1;
            }
        }

        while c_ref.cdigit > 1 && c_ref.mant[(c_ref.cdigit - 1) as usize] == 0 {
            c_ref.cdigit -= 1;
        }
    }

    *pa = c;
}

pub fn numpowi32x(proot: &mut PNumber, mut power: i32) {
    let mut lret = i32tonum(1, BASEX);

    while power > 0 {
        if (power & 1) != 0 {
            mulnumx(&mut lret, std::rc::Rc::clone(proot));
        }

        let proot_clone = std::rc::Rc::clone(proot);
        mulnumx(proot, proot_clone);

        power >>= 1;
    }

    *proot = lret;
}

pub fn divnumx(pa: &mut PNumber, b: PNumber, precision: i32) {
    let b_is_not_one = {
        let b_ref = b.borrow();
        b_ref.cdigit > 1 || b_ref.mant[0] != 1 || b_ref.exp != 0
    };

    if b_is_not_one {
        let a_is_not_one = {
            let a_ref = pa.borrow();
            a_ref.cdigit > 1 || a_ref.mant[0] != 1 || a_ref.exp != 0
        };

        if a_is_not_one {
            _divnumx(pa, b, precision);
        } else {
            let sign = pa.borrow().sign;
            *pa = dupnum(&b);
            pa.borrow_mut().sign *= sign;
        }
    } else {
        let b_sign = b.borrow().sign;
        pa.borrow_mut().sign *= b_sign;
    }
}

pub fn _divnumx(pa: &mut PNumber, b: PNumber, precision: i32) {
    let a = std::rc::Rc::clone(pa);
    
    let mut thismax = precision + 2; 
    if thismax < a.borrow().cdigit {
        thismax = a.borrow().cdigit;
    }
    if thismax < b.borrow().cdigit {
        thismax = b.borrow().cdigit;
    }

    let c = createnum(thismax + 1);
    c.borrow_mut().exp = (a.borrow().cdigit + a.borrow().exp) - (b.borrow().cdigit + b.borrow().exp) + 1;
    c.borrow_mut().sign = a.borrow().sign * b.borrow().sign;

    let mut ptrc = thismax;
    let mut cdigits = 0;

    let mut rem = dupnum(&a);
    rem.borrow_mut().sign = b.borrow().sign;
    rem.borrow_mut().exp = b.borrow().cdigit + b.borrow().exp - rem.borrow().cdigit;

    while cdigits < thismax && !zernum(&rem) {
        let mut digit = 0;
        c.borrow_mut().mant[ptrc as usize] = 0;
        
        while !lessnum(&rem, &b) {
            digit = 1;
            let mut tmp = dupnum(&b);
            let mut lasttmp = i32tonum(0, BASEX);
            
            while lessnum(&tmp, &rem) {
                lasttmp = dupnum(&tmp);
                let tmp_clone = std::rc::Rc::clone(&tmp);
                addnum(&mut tmp, tmp_clone, BASEX);
                digit *= 2;
            }
            
            if lessnum(&rem, &tmp) {
                digit /= 2;
                tmp = lasttmp;
            }

            tmp.borrow_mut().sign *= -1;
            addnum(&mut rem, tmp, BASEX);
            c.borrow_mut().mant[ptrc as usize] |= digit;
        }
        
        rem.borrow_mut().exp += 1;
        ptrc -= 1;
        cdigits += 1;
    }

    if cdigits == 0 {
        c.borrow_mut().exp = 0;
        c.borrow_mut().cdigit = 1;
    } else {
        {
            let mut c_ref = c.borrow_mut();
            let start = (ptrc + 1) as usize;
            for i in 0..(cdigits as usize) {
                c_ref.mant[i] = c_ref.mant[start + i];
            }
            c_ref.cdigit = cdigits;
            c_ref.exp -= cdigits;
            while c_ref.cdigit > 1 && c_ref.mant[(c_ref.cdigit - 1) as usize] == 0 {
                c_ref.cdigit -= 1;
            }
        }
    }

    *pa = c;
}
