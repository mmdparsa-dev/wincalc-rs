use std::cell::RefCell;
use std::rc::Rc;
use std::cmp::{max, min};
use crate::ratpack::types::*;
use crate::ratpack::*;

const MAX_ZEROS_AFTER_DECIMAL: i32 = 2;
const DIGITS: &[u8] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_@";

thread_local! {
    pub static G_DECIMAL_SEPARATOR: RefCell<char> = RefCell::new('.');
}

pub fn set_decimal_separator(decimal_separator: char) {
    G_DECIMAL_SEPARATOR.with(|s| *s.borrow_mut() = decimal_separator);
}

pub fn _dupnum(dest: &mut PNumber, src: &PNumber) {
    let src_ref = src.borrow();
    *dest = Rc::new(RefCell::new(Number {
        sign: src_ref.sign,
        exp: src_ref.exp,
        mant: src_ref.mant.clone(),
    }));
}

pub fn _destroynum(_pnum: Option<&mut PNumber>) {
    // In Rust, memory management is handled by Rc/RefCell
}

pub fn _destroyrat(_prat: Option<&mut PRat>) {
    // In Rust, memory management is handled by Rc/RefCell
}

pub fn _createnum(size: u32) -> PNumber {
    let mut mant = Vec::with_capacity(size as usize);
    mant.resize(size as usize, 0);
    Rc::new(RefCell::new(Number {
        sign: 1,
        exp: 0,
        mant,
    }))
}

pub fn _createrat() -> PRat {
    Rat::new()
}

pub fn numtorat(pin: &PNumber, radix: u32) -> PRat {
    let mut pn_radixn = Number::new();
    _dupnum(&mut pn_radixn, pin);

    let mut qn_radixn = i32tonum(1, radix);

    // Ensure p and q start out as integers.
    {
        let mut pn_radixn_ref = pn_radixn.borrow_mut();
        if pn_radixn_ref.exp < 0 {
            qn_radixn.borrow_mut().exp -= pn_radixn_ref.exp;
            pn_radixn_ref.exp = 0;
        }
    }

    let pout = _createrat();

    {
        let mut pout_ref = pout.borrow_mut();
        pout_ref.pp = numtonRadixx(&pn_radixn, radix);
        pout_ref.pq = numtonRadixx(&qn_radixn, radix);
    }

    pout
}

pub fn nRadixxtonum(a: &PNumber, radix: u32, precision: i32) -> PNumber {
    let mut sum = i32tonum(0, radix);
    let mut powofn_radix = Ui32tonum(BASEX, radix);

    let a_ref = a.borrow();
    let mut cdigits = precision + 1;
    if cdigits > a_ref.mant.len() as i32 {
        cdigits = a_ref.mant.len() as i32;
    }

    numpowi32(&mut powofn_radix, a_ref.exp + (a_ref.mant.len() as i32 - cdigits), radix, precision);

    let mut index = a_ref.mant.len() as i32 - 1;
    while cdigits > 0 {
        let mut bitmask = BASEX / 2;
        while bitmask > 0 {
            addnum(&mut sum, &sum, radix);
            if (a_ref.mant[index as usize] & bitmask) != 0 {
                sum.borrow_mut().mant[0] |= 1;
            }
            bitmask /= 2;
        }
        index -= 1;
        cdigits -= 1;
    }

    mulnum(&mut sum, &powofn_radix, radix);

    sum.borrow_mut().sign = a_ref.sign;
    sum
}

pub fn numtonRadixx(a: &PNumber, radix: u32) -> PNumber {
    let mut pnumret = i32tonum(0, BASEX);
    let mut num_radix = i32tonum(radix as i32, BASEX);
    
    let a_ref = a.borrow();
    let mut index = a_ref.mant.len() as i32 - 1;

    for _ in 0..a_ref.mant.len() {
        mulnumx(&mut pnumret, &num_radix);
        let thisdigit = Ui32tonum(a_ref.mant[index as usize], BASEX);
        addnum(&mut pnumret, &thisdigit, BASEX);
        index -= 1;
    }

    numpowi32x(&mut num_radix, a_ref.exp);
    mulnumx(&mut pnumret, &num_radix);

    pnumret.borrow_mut().sign = a_ref.sign;
    pnumret
}

pub fn string_to_rat(mantissa_is_negative: bool, mantissa: &str, exponent_is_negative: bool, exponent: &str, radix: u32, precision: i32) -> Option<PRat> {
    let mut result_rat: PRat;

    if mantissa.is_empty() {
        if exponent.is_empty() {
            result_rat = Rat::new();
            _dupnum(&mut result_rat.borrow_mut().pp, &rat_zero().borrow().pp);
            _dupnum(&mut result_rat.borrow_mut().pq, &rat_zero().borrow().pq);
        } else {
            result_rat = Rat::new();
            _dupnum(&mut result_rat.borrow_mut().pp, &rat_one().borrow().pp);
            _dupnum(&mut result_rat.borrow_mut().pq, &rat_one().borrow().pq);
        }
    } else {
        let pnummant = string_to_number(mantissa, radix, precision)?;
        result_rat = numtorat(&pnummant, radix);
    }

    let mut expt = 0;
    if !exponent.is_empty() {
        let num_exp = string_to_number(exponent, radix, precision)?;
        expt = numtoi32(&num_exp, radix);
    }

    let mut pnumexp = i32tonum(radix as i32, BASEX);
    numpowi32x(&mut pnumexp, expt.abs());

    let pratexp = _createrat();
    _dupnum(&mut pratexp.borrow_mut().pp, &pnumexp);
    pratexp.borrow_mut().pq = i32tonum(1, BASEX);

    if exponent_is_negative {
        divrat(&mut result_rat, &pratexp, precision);
    } else if expt > 0 {
        mulrat(&mut result_rat, &pratexp, precision);
    }

    if mantissa_is_negative {
        result_rat.borrow_mut().pp.borrow_mut().sign *= -1;
    }

    Some(result_rat)
}

const DP: usize = 0;
const ZR: usize = 1;
const NZ: usize = 2;
const SG: usize = 3;
const EX: usize = 4;

const START: usize = 0;
const MANTS: usize = 1;
const LZ: usize = 2;
const LZDP: usize = 3;
const LD: usize = 4;
const DZ: usize = 5;
const DD: usize = 6;
const DDP: usize = 7;
const EXPB: usize = 8;
const EXPS: usize = 9;
const EXPD: usize = 10;
const EXPBZ: usize = 11;
const EXPSZ: usize = 12;
const EXPDZ: usize = 13;
const ERR: usize = 14;

const MACHINE: [[usize; EX + 1]; ERR + 1] = [
    [LZDP, LZ, LD, MANTS, ERR],
    [LZDP, LZ, LD, ERR, ERR],
    [LZDP, LZ, LD, ERR, EXPBZ],
    [ERR, DZ, DD, ERR, EXPB],
    [DDP, LD, LD, ERR, EXPB],
    [ERR, DZ, DD, ERR, EXPBZ],
    [ERR, DD, DD, ERR, EXPB],
    [ERR, DD, DD, ERR, EXPB],
    [ERR, EXPD, EXPD, EXPS, ERR],
    [ERR, EXPD, EXPD, ERR, ERR],
    [ERR, EXPD, EXPD, ERR, ERR],
    [ERR, EXPDZ, EXPDZ, EXPSZ, ERR],
    [ERR, EXPDZ, EXPDZ, ERR, ERR],
    [ERR, EXPDZ, EXPDZ, ERR, ERR],
    [ERR, ERR, ERR, ERR, ERR],
];

pub fn normalize_char_digit(c: char, radix: u32) -> char {
    if radix as usize >= DIGITS.iter().position(|&x| x == b'A').unwrap() && 
       radix as usize <= DIGITS.iter().position(|&x| x == b'Z').unwrap() {
        return c.to_ascii_uppercase();
    }
    c
}

pub fn string_to_number(number_string: &str, radix: u32, precision: i32) -> Option<PNumber> {
    let mut exp_sign = 1;
    let mut exp_value = 0;

    let pnumret = _createnum(number_string.len() as u32);
    {
        let mut pnum_ref = pnumret.borrow_mut();
        pnum_ref.sign = 1;
        pnum_ref.exp = 0;
        pnum_ref.mant.clear(); // We'll push instead of walking backwards like C++
    }
    
    let mut mant_temp: Vec<u32> = Vec::new();

    let mut state = START;
    let dec_sep = G_DECIMAL_SEPARATOR.with(|s| *s.borrow());

    for c in number_string.chars() {
        let mut cur_char = if c == dec_sep { '.' } else { c };

        let input_type = match cur_char {
            '-' | '+' => SG,
            '.' => DP,
            '0' => ZR,
            '^' | 'e' => {
                if cur_char == '^' || radix == 10 {
                    EX
                } else {
                    NZ
                }
            }
            _ => NZ,
        };

        state = MACHINE[state][input_type];

        match state {
            MANTS => {
                pnumret.borrow_mut().sign = if cur_char == '-' { -1 } else { 1 };
            }
            EXPSZ | EXPS => {
                exp_sign = if cur_char == '-' { -1 } else { 1 };
            }
            EXPDZ | EXPD => {
                cur_char = normalize_char_digit(cur_char, radix);
                if let Some(pos) = DIGITS.iter().position(|&x| x == cur_char as u8) {
                    exp_value *= radix as i32;
                    exp_value += pos as i32;
                } else {
                    state = ERR;
                }
            }
            LD => {
                pnumret.borrow_mut().exp += 1;
                cur_char = normalize_char_digit(cur_char, radix);
                if let Some(pos) = DIGITS.iter().position(|&x| x == cur_char as u8) {
                    if pos < radix as usize {
                        mant_temp.push(pos as u32);
                        let mut pnum_ref = pnumret.borrow_mut();
                        pnum_ref.exp -= 1;
                    } else {
                        state = ERR;
                    }
                } else {
                    state = ERR;
                }
            }
            DD => {
                cur_char = normalize_char_digit(cur_char, radix);
                if let Some(pos) = DIGITS.iter().position(|&x| x == cur_char as u8) {
                    if pos < radix as usize {
                        mant_temp.push(pos as u32);
                        let mut pnum_ref = pnumret.borrow_mut();
                        pnum_ref.exp -= 1;
                    } else {
                        state = ERR;
                    }
                } else {
                    state = ERR;
                }
            }
            DZ => {
                pnumret.borrow_mut().exp -= 1;
            }
            LZ | LZDP | DDP => {}
            _ => {}
        }
    }

    mant_temp.reverse();
    pnumret.borrow_mut().mant = mant_temp;

    if state == DZ || state == EXPDZ {
        let mut pnum_ref = pnumret.borrow_mut();
        pnum_ref.mant = vec![0];
        pnum_ref.exp = 0;
        pnum_ref.sign = 1;
    } else {
        let len = pnumret.borrow().mant.len();
        while pnumret.borrow().mant.len() < number_string.len() {
            pnumret.borrow_mut().mant.insert(0, 0);
            pnumret.borrow_mut().exp -= 1;
        }

        pnumret.borrow_mut().exp += exp_sign * exp_value;
    }

    if pnumret.borrow().mant.is_empty() {
        return None;
    }

    stripzeroesnum(&pnumret, precision);

    Some(pnumret)
}

pub fn i32torat(ini32: i32) -> PRat {
    let pratret = _createrat();
    pratret.borrow_mut().pp = i32tonum(ini32, BASEX);
    pratret.borrow_mut().pq = i32tonum(1, BASEX);
    pratret
}

pub fn Ui32torat(inui32: u32) -> PRat {
    let pratret = _createrat();
    pratret.borrow_mut().pp = Ui32tonum(inui32, BASEX);
    pratret.borrow_mut().pq = i32tonum(1, BASEX);
    pratret
}

pub fn i32tonum(mut ini32: i32, radix: u32) -> PNumber {
    let pnumret = Number::new();
    let mut mant = Vec::new();
    
    let sign = if ini32 < 0 {
        ini32 = -ini32;
        -1
    } else {
        1
    };
    
    let mut value = ini32 as u32;
    loop {
        mant.push(value % radix);
        value /= radix;
        if value == 0 {
            break;
        }
    }
    
    {
        let mut pnum_ref = pnumret.borrow_mut();
        pnum_ref.mant = mant;
        pnum_ref.sign = sign;
        pnum_ref.exp = 0;
    }
    
    pnumret
}

pub fn Ui32tonum(mut inui32: u32, radix: u32) -> PNumber {
    let pnumret = Number::new();
    let mut mant = Vec::new();
    
    loop {
        mant.push(inui32 % radix);
        inui32 /= radix;
        if inui32 == 0 {
            break;
        }
    }
    
    {
        let mut pnum_ref = pnumret.borrow_mut();
        pnum_ref.mant = mant;
        pnum_ref.sign = 1;
        pnum_ref.exp = 0;
    }
    
    pnumret
}

pub fn rattoi32(prat: &PRat, radix: u32, precision: i32) -> Result<i32, String> {
    if rat_gt(prat, &rat_max_i32(), precision) || rat_lt(prat, &rat_min_i32(), precision) {
        return Err("CALC_E_DOMAIN".to_string());
    }

    let pint = _createrat();
    _dupnum(&mut pint.borrow_mut().pp, &prat.borrow().pp);
    _dupnum(&mut pint.borrow_mut().pq, &prat.borrow().pq);

    intrat(&mut pint.clone(), radix, precision);
    
    let pp = pint.borrow().pp.clone();
    let pq = pint.borrow().pq.clone();
    let mut pp_clone = pp.clone();
    divnumx(&mut pp_clone, &pq, precision);
    pint.borrow_mut().pp = pp_clone;
    
    pint.borrow_mut().pq = num_one();

    Ok(numtoi32(&pint.borrow().pp, BASEX))
}

pub fn rattoUi32(prat: &PRat, radix: u32, precision: i32) -> Result<u32, String> {
    if rat_gt(prat, &rat_dword(), precision) || rat_lt(prat, &rat_zero(), precision) {
        return Err("CALC_E_DOMAIN".to_string());
    }

    let pint = _createrat();
    _dupnum(&mut pint.borrow_mut().pp, &prat.borrow().pp);
    _dupnum(&mut pint.borrow_mut().pq, &prat.borrow().pq);

    intrat(&mut pint.clone(), radix, precision);
    
    let pp = pint.borrow().pp.clone();
    let pq = pint.borrow().pq.clone();
    let mut pp_clone = pp.clone();
    divnumx(&mut pp_clone, &pq, precision);
    pint.borrow_mut().pp = pp_clone;
    
    pint.borrow_mut().pq = num_one();

    Ok(numtoi32(&pint.borrow().pp, BASEX) as u32)
}

pub fn rattoUi64(prat: &PRat, radix: u32, precision: i32) -> Result<u64, String> {
    let mut pint = _createrat();
    _dupnum(&mut pint.borrow_mut().pp, &prat.borrow().pp);
    _dupnum(&mut pint.borrow_mut().pq, &prat.borrow().pq);
    
    andrat(&mut pint, &rat_dword(), radix, precision);
    let lo = rattoUi32(&pint, radix, precision)?;

    let mut pint = _createrat();
    _dupnum(&mut pint.borrow_mut().pp, &prat.borrow().pp);
    _dupnum(&mut pint.borrow_mut().pq, &prat.borrow().pq);
    
    let prat32 = i32torat(32);
    rshrat(&mut pint, &prat32, radix, precision);
    intrat(&mut pint, radix, precision);
    andrat(&mut pint, &rat_dword(), radix, precision);
    let hi = rattoUi32(&pint, radix, precision)?;

    Ok(((hi as u64) << 32) | (lo as u64))
}

pub fn numtoi32(pnum: &PNumber, radix: u32) -> i32 {
    let mut lret = 0;
    let pnum_ref = pnum.borrow();
    
    let mut expt = pnum_ref.exp;
    let mut index = pnum_ref.mant.len() as i32 - 1;
    let mut length = pnum_ref.mant.len() as i32;

    while length > 0 && length + expt > 0 {
        lret *= radix as i32;
        lret += pnum_ref.mant[index as usize] as i32;
        index -= 1;
        length -= 1;
    }

    while expt > 0 {
        lret *= radix as i32;
        expt -= 1;
    }
    
    lret *= pnum_ref.sign;
    lret
}

pub fn stripzeroesnum(pnum: &PNumber, starting: i32) -> bool {
    let mut fstrip = false;
    let mut pnum_ref = pnum.borrow_mut();
    
    let mut cdigits = pnum_ref.mant.len() as i32;
    let mut offset = 0;
    
    if cdigits > starting {
        offset = cdigits - starting;
        cdigits = starting;
    }
    
    while cdigits > 0 && pnum_ref.mant[offset as usize] == 0 {
        offset += 1;
        cdigits -= 1;
        fstrip = true;
    }
    
    if fstrip {
        pnum_ref.mant.drain(0..offset as usize);
        pnum_ref.exp += offset;
    }
    
    fstrip
}

pub fn NumberToString(pnum: &mut PNumber, mut format: NumberFormat, radix: u32, precision: i32) -> String {
    stripzeroesnum(pnum, precision + 2);
    let length = pnum.borrow().mant.len() as i32;
    let mut exponent = pnum.borrow().exp + length;

    let old_format = format;
    if exponent > precision && format == NumberFormat::Float {
        format = NumberFormat::Scientific;
    }

    let mut calc_length = length;
    if calc_length > precision {
        calc_length = precision;
    }

    let mut round: Option<PNumber> = None;
    if !zernum(pnum) && (pnum.borrow().mant.len() as i32 >= precision || (calc_length - exponent > precision && exponent >= -MAX_ZEROS_AFTER_DECIMAL)) {
        let mut r = i32tonum(radix as i32, radix);
        divnum(&mut r, &num_two(), radix, precision);

        if exponent > 0 || format == NumberFormat::Float {
            r.borrow_mut().exp = pnum.borrow().exp + pnum.borrow().mant.len() as i32 - r.borrow().mant.len() as i32 - precision;
        } else {
            r.borrow_mut().exp = pnum.borrow().exp + pnum.borrow().mant.len() as i32 - r.borrow().mant.len() as i32 - precision - exponent;
            calc_length = precision + exponent;
        }
        r.borrow_mut().sign = pnum.borrow().sign;
        round = Some(r);
    }

    if format == NumberFormat::Float {
        if (calc_length - exponent > precision) || (exponent > precision + 3) {
            if exponent >= -MAX_ZEROS_AFTER_DECIMAL {
                if let Some(ref r) = round {
                    r.borrow_mut().exp -= exponent;
                }
                calc_length = precision + exponent;
            } else {
                format = NumberFormat::Scientific;
            }
        } else if calc_length + exponent.abs() < precision {
            if let Some(ref r) = round {
                r.borrow_mut().exp -= exponent;
            }
        }
    }

    if let Some(r) = round {
        addnum(pnum, &r, radix);
        let offset = (pnum.borrow().mant.len() as i32 + pnum.borrow().exp) - (r.borrow().mant.len() as i32 + r.borrow().exp);
        if stripzeroesnum(pnum, offset) {
            return NumberToString(pnum, old_format, radix, precision);
        }
    } else {
        stripzeroesnum(pnum, precision);
    }

    let mut use_sci_form = false;
    let mut eout = exponent - 1;
    
    if format == NumberFormat::Scientific || format == NumberFormat::Engineering {
        use_sci_form = true;
        if eout != 0 {
            if format == NumberFormat::Engineering {
                exponent = eout % 3;
                eout -= exponent;
                exponent += 1;
                
                if exponent < 0 {
                    exponent += 3;
                    eout -= 3;
                }
            } else {
                exponent = 1;
            }
        }
    } else {
        eout = 0;
    }

    let mut result = String::new();

    if pnum.borrow().sign == -1 && calc_length > 0 {
        result.push('-');
    }

    if exponent <= 0 && !use_sci_form {
        result.push('0');
        result.push(G_DECIMAL_SEPARATOR.with(|s| *s.borrow()));
    }

    while exponent < 0 {
        result.push('0');
        exponent += 1;
    }

    let mut pmant_idx = pnum.borrow().mant.len() as i32 - 1;
    
    while calc_length > 0 {
        exponent -= 1;
        if pmant_idx >= 0 {
            result.push(DIGITS[pnum.borrow().mant[pmant_idx as usize] as usize] as char);
            pmant_idx -= 1;
        } else {
            result.push('0');
        }
        calc_length -= 1;

        if exponent == 0 {
            result.push(G_DECIMAL_SEPARATOR.with(|s| *s.borrow()));
        }
    }

    while exponent > 0 {
        result.push('0');
        exponent -= 1;
        if exponent == 0 {
            result.push(G_DECIMAL_SEPARATOR.with(|s| *s.borrow()));
        }
    }

    if use_sci_form {
        result.push(if radix == 10 { 'e' } else { '^' });
        result.push(if eout < 0 { '-' } else { '+' });
        eout = eout.abs();
        
        let mut exp_string = String::new();
        loop {
            exp_string.push(DIGITS[(eout % radix as i32) as usize] as char);
            eout /= radix as i32;
            if eout == 0 {
                break;
            }
        }
        result.extend(exp_string.chars().rev());
    }

    if result.ends_with(G_DECIMAL_SEPARATOR.with(|s| *s.borrow())) {
        result.pop();
    }

    result
}

pub fn RatToString(prat: &mut PRat, format: NumberFormat, radix: u32, precision: i32) -> String {
    let mut p = RatToNumber(prat, radix, precision);
    NumberToString(&mut p, format, radix, precision)
}

pub fn RatToNumber(prat: &PRat, radix: u32, precision: i32) -> PNumber {
    let temprat = _createrat();
    _dupnum(&mut temprat.borrow_mut().pp, &prat.borrow().pp);
    _dupnum(&mut temprat.borrow_mut().pq, &prat.borrow().pq);
    
    let scaleby = max(min(temprat.borrow().pp.borrow().exp, temprat.borrow().pq.borrow().exp), 0);
    
    temprat.borrow_mut().pp.borrow_mut().exp -= scaleby;
    temprat.borrow_mut().pq.borrow_mut().exp -= scaleby;

    let mut p = nRadixxtonum(&temprat.borrow().pp, radix, precision);
    let q = nRadixxtonum(&temprat.borrow().pq, radix, precision);

    divnum(&mut p, &q, radix, precision);

    p
}

pub fn flatrat(prat: &mut PRat, radix: u32, precision: i32) {
    let pnum = RatToNumber(prat, radix, precision);
    *prat = numtorat(&pnum, radix);
}

pub fn gcd(a: &PNumber, b: &PNumber) -> PNumber {
    if zernum(a) {
        let mut ret = Number::new();
        _dupnum(&mut ret, b);
        return ret;
    } else if zernum(b) {
        let mut ret = Number::new();
        _dupnum(&mut ret, a);
        return ret;
    }

    let mut larger = Number::new();
    let mut smaller = Number::new();

    if lessnum(a, b) {
        _dupnum(&mut larger, b);
        _dupnum(&mut smaller, a);
    } else {
        _dupnum(&mut larger, a);
        _dupnum(&mut smaller, b);
    }

    while !zernum(&smaller) {
        remnum(&mut larger, &smaller, BASEX);
        let r = larger;
        larger = smaller;
        smaller = r;
    }
    larger
}

pub fn i32factnum(mut ini32: i32, radix: u32) -> PNumber {
    let mut lret = i32tonum(1, radix);

    while ini32 > 0 {
        let tmp = i32tonum(ini32, radix);
        ini32 -= 1;
        mulnum(&mut lret, &tmp, radix);
    }
    lret
}

pub fn i32prodnum(mut start: i32, stop: i32, radix: u32) -> PNumber {
    let mut lret = i32tonum(1, radix);

    while start <= stop {
        if start != 0 {
            let tmp = i32tonum(start, radix);
            mulnum(&mut lret, &tmp, radix);
        }
        start += 1;
    }
    lret
}

pub fn numpowi32(proot: &mut PNumber, mut power: i32, radix: u32, precision: i32) {
    let mut lret = i32tonum(1, radix);

    while power > 0 {
        if (power & 1) != 0 {
            mulnum(&mut lret, proot, radix);
        }
        
        let mut temp_root = Number::new();
        _dupnum(&mut temp_root, proot);
        mulnum(proot, &temp_root, radix);
        
        TRIMNUM(proot, precision);
        power >>= 1;
    }
    
    *proot = lret;
}

pub fn ratpowi32(proot: &mut PRat, mut power: i32, precision: i32) {
    if power < 0 {
        ratpowi32(proot, -power, precision);
        let temp = proot.borrow().pp.clone();
        proot.borrow_mut().pp = proot.borrow().pq.clone();
        proot.borrow_mut().pq = temp;
    } else {
        let mut lret = i32torat(1);

        while power > 0 {
            if (power & 1) != 0 {
                let mut lret_pp = lret.borrow().pp.clone();
                mulnumx(&mut lret_pp, &proot.borrow().pp);
                lret.borrow_mut().pp = lret_pp;
                
                let mut lret_pq = lret.borrow().pq.clone();
                mulnumx(&mut lret_pq, &proot.borrow().pq);
                lret.borrow_mut().pq = lret_pq;
            }
            
            let temp_root = _createrat();
            _dupnum(&mut temp_root.borrow_mut().pp, &proot.borrow().pp);
            _dupnum(&mut temp_root.borrow_mut().pq, &proot.borrow().pq);
            
            mulrat(proot, &temp_root, precision);
            trimit(&mut lret, precision);
            trimit(proot, precision);
            power >>= 1;
        }
        *proot = lret;
    }
}
