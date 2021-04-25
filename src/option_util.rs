use std::f64::consts::{PI, E, SQRT_2, FRAC_1_SQRT_2};
use libm::erfc;

fn get_d1(s: f64, k: f64, t: f64, b: f64, vol: f64) -> f64 {
    return ((s / k).log(E) + (b + vol.powi(2) / 2.0) * t) / (vol.max(0.0000000000001) * t.sqrt());
}

fn get_d2(t: f64, vol: f64, d1: f64) -> f64 {
    return d1 - vol * t.sqrt();
}

/// 标准正太分布概率密度
fn nd(value: f64) -> f64 {
    return (1.0 / (2.0 * PI).sqrt()) * (-value.powi(2) / 2.0).exp();
}

/// 正太分布
fn normal_cdf(value: f64) -> f64 {
    return 0.5 * erfc(-value * FRAC_1_SQRT_2);
}

/// Black-Scholes期权定价模型
pub fn black_scholes_calc(s: f64, k: f64, t: f64, r: f64, b: f64, vol: f64, call_put: char) -> f64 {
    let mut bs = 0.0;

    let d1 = get_d1(s, k, t, b, vol);
    let d2 = get_d2(t, vol, d1);

    if call_put == 'c' {
        let a = ((b - r) * t).exp() * normal_cdf(d1);
        let b = (-r * t).exp() * normal_cdf(d2);
        bs = s * a - k * b;
    } else if call_put == 'p' {
        bs = k * (-r * t).exp() * normal_cdf(-d2) - s * ((b - r) * t).exp() * normal_cdf(-d1);
    }

    return bs;
}

/// 隐含波动率
///
/// # Arguments
///
/// * `s` - 标的资产价格
/// * `k` - 行权价
/// * `t` - 到期日
/// * `r` - 无风险利率
/// * `b` - 股票红利率
/// * `target` -
/// * `call_put` - 认沽/认购
///
pub fn implied_volatility(s: f64, k: f64, t: f64, r: f64, b: f64, target: f64, call_put: char) -> f64 {
    let mut high = 5.0;
    let mut low = 0.0;

    while high - low > 0.00000001 {
        let vol = (high + low) / 2.0;
        let bs = black_scholes_calc(s, k, t, r, b, vol, call_put);

        if bs > target {
            high = vol;
        } else {
            low = vol;
        }
    }

    return (high + low) / 2.0;
}

pub fn delta(s: f64, k: f64, t: f64, r: f64, b: f64, vol: f64, cp: char) -> f64 {
    let mut delta = 0.0;
    let d1 = get_d1(s, k, t, b, vol);

    if cp == 'c' {
        delta = ((b - r) * t).exp() * normal_cdf(d1);
    } else if cp == 'p' {
        delta = ((b - r) * t).exp() * (normal_cdf(d1) - 1.0);
    }

    return delta;
}

pub fn gamma(s: f64, k: f64, t: f64, r: f64, b: f64, vol: f64) -> f64 {
    let d1 = get_d1(s, k, t, b, vol);
    return ((b - r) * t).exp() * nd(d1) / (s * vol.max(0.0000000000001) * t.sqrt());
}

pub fn vega(s: f64, k: f64, t: f64, r: f64, b: f64, vol: f64) -> f64 {
    let d1 = get_d1(s, k, t, b, vol);
    return s * ((b - r) * t).exp() * t.sqrt() * nd(d1);
}

pub fn theta(s: f64, k: f64, t: f64, r: f64, b: f64, vol: f64, target: f64, cp: char) -> f64 {
    let mut theta = 0.0;
    let d1 = get_d1(s, k, t, b, vol);
    let d2 = get_d2(t, vol, d1);

    if cp == 'c' {
        theta = -s * ((b - r) * t).exp() * nd(d1) * vol / (2.0 * t.sqrt()) - (b - r) * s
            * ((b - r) * t).exp() * normal_cdf(d1) - r * k * (-r * t).exp()
            * normal_cdf(d2);
    } else if cp == 'p' {
        theta = -s * (-b * t).exp() * nd(d1) * vol / (2.0 * t.sqrt()) - (b - r) * s
            * ((b - r) * t).exp() * normal_cdf(-d1) - r * k * (-r * t).exp()
            * normal_cdf(d2);
    }

    if (theta * t + target) > 0.0 {
        return theta;
    }

    return (black_scholes_calc(s, k, t - (1.0 / 246.0), r, b, vol, cp)
        - black_scholes_calc(s, k, t, r, b, vol, cp)) * 246.0;
}

pub fn rho(s: f64, k: f64, t: f64, r: f64, b: f64, vol: f64, cp: char) -> f64 {
    let mut rho = 0.0;
    let d1 = get_d1(s, k, t, b, vol);
    let d2 = get_d2(t, vol, d1);

    if cp == 'c' {
        rho = k * t * (-r * t).exp() * normal_cdf(d2);
    } else if cp == 'p' {
        rho = -k * t * (-r * t).exp() * normal_cdf(-d2);
    }

    return rho;
}

pub fn effective_time_etf(t: f64) -> f64 {
    let mut _t = t;
    let mut ret = 0.0;
    _t /= 1000.0;

    let t1 = 9.0 * 3.6 + 30.0 * 0.06;
    let t2 = 11.0 * 3.6 + 30.0 * 0.06;
    let t3 = 13.0 * 3.6;
    let t4 = 15.0 * 3.6;

    let total_time = 4.0 * 3.6;

    if _t >= t1 && _t < t3 {
        ret = 1.0 - (_t - t1).min(t2 - t1) / total_time;
    } else if _t >= t3 && _t < t4 {
        ret = 1.0 - (_t - t3 + 2.0 * 3.6) / total_time;
    } else if _t < t1 {
        ret = 1.0;
    }

    return ret;
}

pub struct OptionUtil {
    s: f64,
    k: f64,
    t: f64,
    r: f64,
    b: f64,
    vol: f64,
    target: f64,
    call_put: char,
}

impl OptionUtil {
    pub fn new(spot_price: f64, strike_price: f64, maturity_time: f64, risk_free_rate: f64,
               dividend: f64, vol: f64, target: f64, call_put: char) -> OptionUtil {
        OptionUtil {
            s: spot_price,
            k: strike_price,
            t: maturity_time,
            r: risk_free_rate,
            b: dividend,
            vol,
            target,
            call_put,
        }
    }

    pub fn delta(&self) -> f64 {
        return delta(self.s, self.k, self.t, self.r, self.b, self.vol, self.call_put);
    }

    pub fn iv(&self) -> f64 {
        return implied_volatility(self.s, self.k, self.t, self.r, self.b, self.target, self.call_put);
    }

    pub fn gamma(&self) -> f64 {
        return gamma(self.s, self.k, self.t, self.r, self.b, self.vol);
    }

    pub fn vega(&self) -> f64 {
        return vega(self.s, self.k, self.t, self.r, self.b, self.vol);
    }

    pub fn theta(&self) -> f64 {
        return theta(self.s, self.k, self.t, self.r, self.b, self.vol, self.target, self.call_put);
    }

    pub fn rho(&self) -> f64 {
        return rho(self.s, self.k, self.t, self.r, self.b, self.vol, self.call_put);
    }
}
