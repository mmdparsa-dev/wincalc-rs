use crate::ratpack::*;
use std::cmp::{max, min};
use std::collections::VecDeque;
use std::rc::Rc;

pub fn addnum(pa: &mut PNumber, b: PNumber, radix: u32) {
    let b_is_nonzero = {
        let b_ref = b.borrow();
        b_ref.cdigit > 1 || b_ref.mant[0] != 0
    };
    if b_is_nonzero {
        let a_is_nonzero = {
            let a_ref = pa.borrow();
            a_ref.cdigit > 1 || a_ref.mant[0] != 0
        };
        if a_is_nonzero {
            _addnum(pa, b, radix);
        } else {
            *pa = dupnum(&b);
        }
    }
}

fn _addnum(pa: &mut PNumber, b: PNumber, radix: u32) {
    let a = Rc::clone(pa);

    let (a_cdigit, a_exp, a_sign) = {
        let a_ref = a.borrow();
        (a_ref.cdigit, a_ref.exp, a_ref.sign)
    };
    let (b_cdigit, b_exp, b_sign) = {
        let b_ref = b.borrow();
        (b_ref.cdigit, b_ref.exp, b_ref.sign)
    };

    let mut cdigits = max(a_cdigit + a_exp, b_cdigit + b_exp) - min(a_exp, b_exp);

    let c = createnum(cdigits + 1);
    let c_exp = min(a_exp, b_exp);
    let c_cdigit = cdigits;

    let mut mexp = c_exp;
    let mut cy: u32 = 0;
    let mut fcompla = false;
    let mut fcomplb = false;

    if a_sign != b_sign {
        cy = 1;
        fcompla = a_sign == -1;
        fcomplb = b_sign == -1;
    }

    let mut a_idx = 0;
    let mut b_idx = 0;
    let mut c_idx = 0;

    {
        let a_ref = a.borrow();
        let b_ref = b.borrow();
        let mut c_ref = c.borrow_mut();

        c_ref.exp = c_exp;
        c_ref.cdigit = c_cdigit;

        while cdigits > 0 {
            let mut da: u32 = if mexp >= a_ref.exp && (cdigits + a_ref.exp - c_ref.exp > (c_ref.cdigit - a_ref.cdigit)) {
                let val = a_ref.mant[a_idx as usize];
                a_idx += 1;
                val
            } else {
                0
            };

            let mut db: u32 = if mexp >= b_ref.exp && (cdigits + b_ref.exp - c_ref.exp > (c_ref.cdigit - b_ref.cdigit)) {
                let val = b_ref.mant[b_idx as usize];
                b_idx += 1;
                val
            } else {
                0
            };

            if fcompla {
                da = radix - 1 - da;
            }
            if fcomplb {
                db = radix - 1 - db;
            }

            cy = da + db + cy;
            c_ref.mant[c_idx as usize] = cy % radix;
            c_idx += 1;
            cy /= radix;

            cdigits -= 1;
            mexp += 1;
        }

        if cy != 0 && !(fcompla || fcomplb) {
            c_ref.mant[c_idx as usize] = cy;
            c_idx += 1;
            c_ref.cdigit += 1;
        }

        if !(fcompla || fcomplb) {
            c_ref.sign = a_sign;
        } else {
            if cy != 0 {
                c_ref.sign = 1;
            } else {
                c_ref.sign = -1;
                cy = 1;
                c_idx = 0;
                let mut temp_cdigits = c_ref.cdigit;
                while temp_cdigits > 0 {
                    cy = radix - 1 - c_ref.mant[c_idx as usize] + cy;
                    c_ref.mant[c_idx as usize] = cy % radix;
                    c_idx += 1;
                    cy /= radix;
                    temp_cdigits -= 1;
                }
            }
        }

        while c_ref.cdigit > 1 && c_ref.mant[(c_ref.cdigit - 1) as usize] == 0 {
            c_ref.cdigit -= 1;
        }
    }

    *pa = c;
}

pub fn mulnum(pa: &mut PNumber, b: PNumber, radix: u32) {
    let b_is_one = {
        let b_ref = b.borrow();
        b_ref.cdigit <= 1 && b_ref.mant[0] == 1 && b_ref.exp == 0
    };
    if !b_is_one {
        let a_is_one = {
            let a_ref = pa.borrow();
            a_ref.cdigit <= 1 && a_ref.mant[0] == 1 && a_ref.exp == 0
        };
        if !a_is_one {
            _mulnum(pa, b, radix);
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

fn _mulnum(pa: &mut PNumber, b: PNumber, radix: u32) {
    let a = Rc::clone(pa);

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
                    cy += (c_ref.mant[icdigit as usize] as u64) + (mcy % (radix as u64));
                    c_ref.mant[icdigit as usize] = (cy % (radix as u64)) as u32;
                    icdigit += 1;
                    mcy /= radix as u64;
                    cy /= radix as u64;
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

fn msd(n: &PNumber) -> u32 {
    let n_ref = n.borrow();
    n_ref.mant[(n_ref.cdigit - 1) as usize]
}

pub fn remnum(pa: &mut PNumber, b: PNumber, radix: u32) {
    let mut tmp: PNumber;
    let mut lasttmp: PNumber;

    while !lessnum(pa, &b) {
        tmp = dupnum(&b);
        if lessnum(&tmp, pa) {
            let pa_cdigit = pa.borrow().cdigit;
            let pa_exp = pa.borrow().exp;
            let tmp_cdigit = tmp.borrow().cdigit;
            tmp.borrow_mut().exp = pa_cdigit + pa_exp - tmp_cdigit;
            if msd(pa) <= msd(&tmp) {
                tmp.borrow_mut().exp -= 1;
            }
        }

        lasttmp = i32tonum(0, radix);

        while lessnum(&tmp, pa) {
            lasttmp = dupnum(&tmp);
            let mut tmp_clone = Rc::clone(&tmp);
            addnum(&mut tmp, tmp_clone, radix);
        }

        if lessnum(pa, &tmp) {
            tmp = lasttmp;
        }

        let pa_sign = pa.borrow().sign;
        tmp.borrow_mut().sign = -1 * pa_sign;
        addnum(pa, tmp, radix);
    }
}

pub fn divnum(pa: &mut PNumber, b: PNumber, radix: u32, precision: i32) {
    let b_is_one = {
        let b_ref = b.borrow();
        b_ref.cdigit <= 1 && b_ref.mant[0] == 1 && b_ref.exp == 0
    };

    if !b_is_one {
        _divnum(pa, b, radix, precision);
    } else {
        let b_sign = b.borrow().sign;
        pa.borrow_mut().sign *= b_sign;
    }
}

fn _divnum(pa: &mut PNumber, b: PNumber, radix: u32, precision: i32) {
    let a = Rc::clone(pa);

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

    let mut rem = dupnum(&a);
    let tmp = dupnum(&b);
    tmp.borrow_mut().sign = a.borrow().sign;
    rem.borrow_mut().exp = b.borrow().cdigit + b.borrow().exp - rem.borrow().cdigit;

    let mut number_list = VecDeque::new();
    number_list.push_back(i32tonum(0, radix));

    for _ in 1..radix {
        let mut new_value = dupnum(number_list.front().unwrap());
        addnum(&mut new_value, Rc::clone(&tmp), radix);
        number_list.push_front(new_value);
    }

    let mut ptrc = thismax;
    let mut cdigits = 0;

    while cdigits < thismax && !zernum(&rem) {
        cdigits += 1;
        let mut digit = radix - 1;
        let mut multiple = None;

        for num in &number_list {
            if !lessnum(&rem, num) {
                multiple = Some(Rc::clone(num));
                break;
            }
            digit -= 1;
            if digit == 0 {
                multiple = Some(Rc::clone(num));
                break;
            }
        }

        if let Some(multiple) = multiple {
            if digit != 0 {
                multiple.borrow_mut().sign *= -1;
                addnum(&mut rem, Rc::clone(&multiple), radix);
                multiple.borrow_mut().sign *= -1;
            }
        }
        rem.borrow_mut().exp += 1;
        
        c.borrow_mut().mant[ptrc as usize] = digit;
        ptrc -= 1;
    }

    if cdigits == 0 {
        c.borrow_mut().cdigit = 1;
        c.borrow_mut().exp = 0;
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

pub fn equnum(a: &PNumber, b: &PNumber) -> bool {
    let (a_cdigit, a_exp) = {
        let a_ref = a.borrow();
        (a_ref.cdigit, a_ref.exp)
    };
    let (b_cdigit, b_exp) = {
        let b_ref = b.borrow();
        (b_ref.cdigit, b_ref.exp)
    };

    let diff = (a_cdigit + a_exp) - (b_cdigit + b_exp);
    if diff != 0 {
        return false;
    }

    let cdigits_max = max(a_cdigit, b_cdigit);
    let ccdigits = cdigits_max;

    let a_ref = a.borrow();
    let b_ref = b.borrow();

    let mut a_idx = a_cdigit - 1;
    let mut b_idx = b_cdigit - 1;

    for cdigits in (1..=cdigits_max).rev() {
        let da = if cdigits > (ccdigits - a_cdigit) {
            let val = a_ref.mant[a_idx as usize];
            a_idx = a_idx.saturating_sub(1);
            val
        } else {
            0
        };

        let db = if cdigits > (ccdigits - b_cdigit) {
            let val = b_ref.mant[b_idx as usize];
            b_idx = b_idx.saturating_sub(1);
            val
        } else {
            0
        };

        if da != db {
            return false;
        }
    }

    true
}

pub fn lessnum(a: &PNumber, b: &PNumber) -> bool {
    let (a_cdigit, a_exp) = {
        let a_ref = a.borrow();
        (a_ref.cdigit, a_ref.exp)
    };
    let (b_cdigit, b_exp) = {
        let b_ref = b.borrow();
        (b_ref.cdigit, b_ref.exp)
    };

    let diff = (a_cdigit + a_exp) - (b_cdigit + b_exp);
    if diff < 0 {
        return true;
    }
    if diff > 0 {
        return false;
    }

    let cdigits_max = max(a_cdigit, b_cdigit);
    let ccdigits = cdigits_max;

    let a_ref = a.borrow();
    let b_ref = b.borrow();

    let mut a_idx = a_cdigit - 1;
    let mut b_idx = b_cdigit - 1;

    for cdigits in (1..=cdigits_max).rev() {
        let da = if cdigits > (ccdigits - a_cdigit) {
            let val = a_ref.mant[a_idx as usize];
            a_idx = a_idx.saturating_sub(1);
            val
        } else {
            0
        };

        let db = if cdigits > (ccdigits - b_cdigit) {
            let val = b_ref.mant[b_idx as usize];
            b_idx = b_idx.saturating_sub(1);
            val
        } else {
            0
        };

        if da != db {
            return da < db;
        }
    }

    false
}

pub fn zernum(a: &PNumber) -> bool {
    let a_ref = a.borrow();
    let mut length = a_ref.cdigit;
    let mut pcha = 0;

    while length > 0 {
        if a_ref.mant[pcha] != 0 {
            return false;
        }
        pcha += 1;
        length -= 1;
    }

    true
}
