extern crate statrs;
extern crate libm;
extern crate mysql;
extern crate serde_json;

mod option_util;
mod option_process;

use crate::option_util::{implied_volatility, black_scholes_calc, theta, effective_time_etf, rho,
                         delta, OptionUtil};
use std::fs;
use mysql::serde::__private::ser::constrain;
use serde_json::{Result, Value};


fn main2() {
    let s = 3.44516;
    let k = 3.1;
    let t = 0.040323788;
    let r = 0.0;
    let b = 0.0;
    let target = 0.34645;
    let _vol = 27.65;
    let e_t = 47721.0;

    let ou = OptionUtil::new(s, k, t, r,
                             b, _vol, target, 'p');
}

fn main() -> Result<()>{
    let contents = fs::read_to_string("/Users/x2h1z/CLionProjects/option-process-rust/data/option_chains_2020.json")
        .expect("Something went wrong reading the file");
    let v: Value = serde_json::from_str(contents.as_str())?;

    for item in v.as_object().iter() {
        for (k, v) in item.iter() {

        }
    }

    Ok(())
}