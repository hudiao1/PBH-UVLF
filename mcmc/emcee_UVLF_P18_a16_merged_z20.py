# 导入计算所需包
from pathlib import Path

import numpy as np
import camb
from camb import model
from scipy.interpolate import PchipInterpolator
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.special import erf
from scipy.interpolate import interp1d, RectBivariateSpline
import emcee
from schwimmbad import MPIPool


# Beta function for the dust extinction
def betaAverage(z, MUV):
    if MUV < -19.5:
        return dbetadMUVinterp(z) * (MUV + 19.5) + betainterp(z)
    return (betainterp(z) + 2.33) * np.exp(
        (dbetadMUVinterp(z) * (MUV + 19.5)) / (betainterp(z) + 2.33)
    ) - 2.33


# Dust extinction parameter
def AUV(z, MUV):
    if z < 2.5 or z > 8:
        return 0.0
    sigmabeta = 0.34
    return max(
        0.0,
        4.54
        + 0.2 * np.log(10) * (2.07**2) * (sigmabeta**2)
        + 2.07 * betaAverage(z, MUV),
    )


def D_p(z, Omega_m, Omega_b):
    r = (Omega_m - Omega_b) / Omega_m
    a = 1 / 4 * ((1 + 24 * r) ** 0.5 - 1)
    d2 = (1 + 3 * r * 3401 / (2 * (1 + z) * a)) ** a
    return d2


# PBH power spectrum
def p_PBH(logM_PBH, logf_PBH, Omega_m, Omega_b, h):
    return (10**logM_PBH * 10**logf_PBH) / (
        (Omega_m - Omega_b) * rho_crit * h**2
    )


# top hat window function
def w_f(k, R):
    x = k * R
    return 3 * (np.sin(x) - x * np.cos(x)) / x**3

k_sigma = np.logspace(-5, 1, 1000)  # Mpc^-1 单位
k_sigma_factor = k_sigma**2 / (2.0 * np.pi**2)
def sigma_PBH(R, z, logM_PBH, logf_PBH, Omega_m, Omega_b, h):
    """
    使用全局变量 sigma_cb_interp_func 来获取 sigma_cb(R, z)
    """
    integrand = (
        p_PBH(logM_PBH, logf_PBH, Omega_m, Omega_b, h)
        * D_p(z, Omega_m, Omega_b) ** 2
        * w_f(k_sigma, R) ** 2
        * k_sigma_factor
    )
    sigma_cb_val = sigma_cb_interp_func(R, z)
    return np.sqrt(simpson(integrand, x=k_sigma) + sigma_cb_val**2)


def lambda_of_z_from_loglambda(
    z,
    loglambda_E68,
    loglambda_E9,
    loglambda_E10,
    loglambda_E11,
    loglambda_E12,
    loglambda_E13p5,
):
    """
    分段定义 Eddington ratio：

    6.0   <= z < 8.5    : loglambda_E68
    8.5   <= z < 9.5    : loglambda_E9
    9.5   <= z < 10.5   : loglambda_E10
    10.5  <= z < 11.5   : loglambda_E11
    11.5  <= z < 13.5   : loglambda_E12
    13.5  <= z <= 1000  : loglambda_E13p5
    """
    z = np.asarray(z)

    loglam = np.zeros_like(z, dtype=float)

    loglam[(z >= 6.0) & (z < 8.5)] = loglambda_E68
    loglam[(z >= 8.5) & (z < 9.5)] = loglambda_E9
    loglam[(z >= 9.5) & (z < 10.5)] = loglambda_E10
    loglam[(z >= 10.5) & (z < 11.5)] = loglambda_E11
    loglam[(z >= 11.5) & (z < 13.5)] = loglambda_E12
    loglam[(z >= 13.5) & (z <= 1000.0)] = loglambda_E13p5

    return 10.0 ** loglam

def camb_Hz_Gyr_inv(camb_results, z):
    """
    CAMB H(z): km/s/Mpc -> Gyr^-1
    """

    Hz_km_s_Mpc = camb_results.hubble_parameter(z)

    return Hz_km_s_Mpc * 1.022712165e-3

def build_MBH_growth_factor_from_camb(
    camb_results,
    loglambda_E68,
    loglambda_E9,
    loglambda_E10,
    loglambda_E11,
    loglambda_E12,
    loglambda_E13p5,
    f_duty=1.0,
    t_Sal=0.045,
    z_merge=20.0,
    z_final=6.0,
    n_points=5000,
):
    """返回 z=20 合并后单个大黑洞的无量纲质量增长因子。"""
    z_bh_grid = np.linspace(z_merge, z_final, n_points)

    lambda_grid = lambda_of_z_from_loglambda(
        z_bh_grid,
        loglambda_E68,
        loglambda_E9,
        loglambda_E10,
        loglambda_E11,
        loglambda_E12,
        loglambda_E13p5,
    )

    H_grid = camb_Hz_Gyr_inv(camb_results, z_bh_grid)

    if np.any(H_grid <= 0.0):
        raise ValueError("H(z) contains non-positive values.")

    integrand = (
        f_duty
        * lambda_grid
        / (t_Sal * (1.0 + z_bh_grid) * H_grid)
    )

    cumulative_integral = cumulative_trapezoid(
        integrand,
        z_bh_grid,
        initial=0.0,
    )

    ln_growth = -cumulative_integral

    if np.any(ln_growth > 700):
        raise OverflowError("Black hole growth exponent too large.")

    growth_factor_grid = np.exp(ln_growth)

    if np.any(~np.isfinite(growth_factor_grid)):
        raise FloatingPointError("M_BH growth factor contains non-finite values.")

    growth_factor_of_z = interp1d(
        z_bh_grid[::-1],
        growth_factor_grid[::-1],
        kind="linear",
        bounds_error=True,
    )

    return growth_factor_of_z, z_bh_grid, growth_factor_grid


def merged_pbh_mass(Mh, logf_PBH):
    """z=20 时晕内所有 PBH 瞬时合并后的总质量。"""
    return float(Mh) * np.power(10.0, float(logf_PBH))

def L1_PBH(z,
    loglambda_E,
    logf_PBH,
    Mh,
    MBH_growth_factor_of_z,):
    M_BH_z = (
        merged_pbh_mass(Mh, logf_PBH)
        * float(MBH_growth_factor_of_z(z))
    )
    if not np.isfinite(M_BH_z):
        raise FloatingPointError("Merged M_BH contains a non-finite value.")
    L_b = 3.28e4 * 3.846e33 * np.power(10, loglambda_E) * M_BH_z
    if not np.isfinite(L_b):
        raise FloatingPointError("Merged black-hole luminosity is non-finite.")
    return (
        L_b
        / (
            1.862 * (float(L_b) / 3.846e43) ** (-0.361)
            + 4.870 * (L_b / 3.846e43) ** (-0.0063)
        )
        * 4.8367e-16
    )


# Comoving radial distance
def rcomoving(z, Omega_m, h):
    return (
        c
        * integrator(
            lambda x: 1 / np.sqrt(Omega_m * np.power(1 + x, 3) + 1.0 - Omega_m),
            0.0,
            z,
        )
        / (100.0 * h)
    )


# Gaussian integrator
def integrator(f, a, b):
    sub = (b - a) / 2.0
    add = (b + a) / 2.0
    if sub == 0:
        return 0.0
    return sub * np.dot(f(sub * points + add), weights)


# Sheth-Mo-Tormen HMF
def HMF(z, M, logM_PBH, logf_PBH, h, Omega_b, Omega_m):
    rhoM = np.power(h, 2) * Omega_m * rho_crit
    deltaM = 100
    sigma = sigma_PBH(
        np.power(3.0 * M / (4.0 * np.pi * rhoM), 1.0 / 3),
        z, logM_PBH, logf_PBH, Omega_m, Omega_b, h,
    )
    dsigmadM = (
        sigma_PBH(
            np.power(3.0 * (M + deltaM) / (4.0 * np.pi * rhoM), 1.0 / 3),
            z, logM_PBH, logf_PBH, Omega_m, Omega_b, h,
        )
        - sigma_PBH(
            np.power(3.0 * (M - deltaM) / (4.0 * np.pi * rhoM), 1.0 / 3),
            z, logM_PBH, logf_PBH, Omega_m, Omega_b, h,
        )
    ) / (2.0 * deltaM)
    return (
        -AST
        * np.sqrt(2.0 * aST / np.pi)
        * (
            1.0
            + np.power(
                np.power(sigma, 2) / (aST * np.power(deltaST, 2)),
                pST,
            )
        )
        * (deltaST / sigma)
        * np.exp(-aST * np.power(deltaST, 2) / (2.0 * np.power(sigma, 2)))
        * (rhoM / (M * sigma))
        * dsigmadM
    )


# Obtain average <MUV> from given halo mass Mh
def MUV_from_Mh(
    growth_interp, growthf_interp, hubble_interp,
    z, Mh, logM_PBH, logf_PBH, h, Omega_b, Omega_m,
    loglambda_E68, loglambda_E9, loglambda_E10, loglambda_E11, loglambda_E12,
    alphastar, betastar, epsilonstar_slope, epsilonstar_icept,
    Mc_slope, Mc_icept, Q, MBH_growth_factor_of_z,
):
    z_float = float(z)

    if z_float in [6.0, 7.0, 8.0]:
        loglambda_E = float(loglambda_E68)
    elif z_float == 9.0:
        loglambda_E = float(loglambda_E9)    
    elif z_float == 10.0:
        loglambda_E = float(loglambda_E10)
    elif z_float == 11.0:
        loglambda_E = float(loglambda_E11)
    elif z_float == 12.5:
        loglambda_E = float(loglambda_E12)
    else:
        raise ValueError(
            f"Redshift z={z} is not in the defined ranges (6, 7, 8, 9, 10, 11, 12.5)."
        )

    Mh = float(Mh)
    logM_PBH = float(logM_PBH)
    epsilonstar = 10 ** (
        epsilonstar_slope * np.log10((1 + z) / (1 + 6)) + epsilonstar_icept
    )
    Mc = 10 ** (Mc_slope * np.log10((1 + z) / (1 + 6)) + Mc_icept)

    rhoM = np.power(h, 2) * Omega_m * rho_crit

    sigma1 = sigma_PBH(
        np.power(3.0 * Mh / Q / (4.0 * np.pi * rhoM), 1.0 / 3),
        z, logM_PBH, logf_PBH, Omega_m, Omega_b, h,
    )
    sigma2 = sigma_PBH(
        np.power(3.0 * Mh / (4.0 * np.pi * rhoM), 1.0 / 3),
        z, logM_PBH, logf_PBH, Omega_m, Omega_b, h,
    )

    functionf = 1.0 / np.sqrt(sigma1**2 - sigma2**2)

    # 使用预计算的插值
    Dz = growth_interp(z)
    fz = growthf_interp(z)
    Hz = hubble_interp(z)  # 单位: 1/Mpc

    dgrowthdz = -fz * Dz / (1.0 + z)

    Mhdot = (
        -(1 + z)
        * Hz
        * invMpctoinvYear
        * 1.686
        * np.sqrt(2.0 / np.pi)
        * Mh
        * functionf
        * dgrowthdz
        / Dz**2
    )

    return (
        -2.5
        * np.log10(
            epsilonstar
            * Mhdot
            / ((Mh / Mc) ** alphastar + (Mh / Mc) ** betastar)
            / kappaUV
            + L1_PBH(z,
    loglambda_E,
    logf_PBH,
    Mh,
    MBH_growth_factor_of_z,)
        )
        + 51.63
    )


# Integrated Gaussian distribution
def first_integrand(MUV, width, MUV_av, sigma_MUV):
    return 0.5 * (
        erf((MUV_av - MUV + width / 2.0) / (sigma_MUV * np.sqrt(2)))
        - erf((MUV_av - MUV - width / 2.0) / (sigma_MUV * np.sqrt(2)))
    )


# ============================================================
# 构建 sigma_cb(R, z) 的二维插值
# ============================================================
def build_sigma_cb_interp(camb_results, zs_for_sigma):
    """
    预计算 sigma_cb(R, z) 的二维插值函数。
    
    CAMB get_sigmaR 的正确用法：
    - z_indices=None 时返回所有预计算 redshift 的结果，形状 (n_z, n_R)
    - z_indices=i 时返回第 i 个 redshift 的结果，形状 (n_R,)
    - redshifts 在 CAMB 内部按**降序**存储
    """
    R_arr = np.logspace(-2, 2.5, 200)  # Mpc 单位

    # CAMB 内部将 redshift 按降序排列
    # zs_for_sigma 传入时是升序的，CAMB 内部会反转
    # get_sigmaR 返回的 z 维度也是按 CAMB 内部降序排列的

    # 一次性获取所有 z 的 sigma(R)，返回 shape = (n_z, n_R)
    sigma_all = camb_results.get_sigmaR(
        R_arr, hubble_units=False,
        var1='delta_tot', var2='delta_tot'
    )
    # sigma_all 的 z 维度按 CAMB 内部顺序（降序）排列
    # 获取 CAMB 内部存储的 redshift 顺序
    # transfer_redshifts 是 CAMB 内部降序排列的 redshift 列表
    z_camb = camb_results.transfer_redshifts  # 降序

    # 反转为升序，方便插值
    z_ascending = np.array(z_camb[::-1])
    sigma_ascending = sigma_all[::-1, :]  # 对应反转 z 维度

    # 构建二维插值: log10(R) vs z -> sigma
    logR_arr = np.log10(R_arr)

    # 使用 RectBivariateSpline 进行二维插值
    # 注意：RectBivariateSpline 要求 x, y 严格递增
    interp_2d = RectBivariateSpline(z_ascending, logR_arr, sigma_ascending, kx=3, ky=3)

    def sigma_cb_func(R, z):
        """返回 sigma_cb(R, z)，R 单位 Mpc"""
        logR = np.log10(R)
        val = interp_2d(z, logR, grid=False)
        return float(val)

    return sigma_cb_func


# ============================================================
# 构建增长因子和 Hubble 的插值
# ============================================================
def build_growth_interps(camb_results, pars):
    """
    构建:
    - D(z): scale-independent growth factor (归一化到 z=0 时 D=1)
    - f(z): growth rate f = d ln D / d ln a
    - H(z): Hubble parameter in 1/Mpc units
    """
    z_arr = np.linspace(0, 14, 500)
    h = pars.H0 / 100.0

    # 获取 sigma8(z=0) 用于归一化
    sigma8_0 = camb_results.get_sigmaR(
        8.0 / h, hubble_units=False,
        var1='delta_tot', var2='delta_tot'
    )
    # sigma8_0 shape = (n_z,)，z=0 是最后一个（CAMB 内部降序）
    # 取 z=0 对应的值
    sigma8_0_val = float(sigma8_0[-1])  # 最后一个是 z=0（降序中最小的 z）

    # 对每个 z 计算 D(z) = sigma_cb(R=8/h, z) / sigma_cb(R=8/h, 0)
    # 使用已经构建好的 sigma_cb_interp_func
    D_arr = np.zeros_like(z_arr)
    R8 = 8.0 / h
    for i, z in enumerate(z_arr):
        sig_z = sigma_cb_interp_func(R8, z)
        D_arr[i] = sig_z / sigma8_0_val if sigma8_0_val > 0 else 0.0

    # f(z) = - (1+z) / D * dD/dz，用数值微分
    f_arr = np.zeros_like(z_arr)
    for i in range(1, len(z_arr) - 1):
        dDdz = (D_arr[i + 1] - D_arr[i - 1]) / (z_arr[i + 1] - z_arr[i - 1])
        f_arr[i] = -(1 + z_arr[i]) / D_arr[i] * dDdz if D_arr[i] > 0 else 0.0
    f_arr[0] = f_arr[1]
    f_arr[-1] = f_arr[-2]

    # H(z) in 1/Mpc units: H_camb [km/s/Mpc] -> H [1/Mpc] = H_camb / c
    H_arr = np.array([camb_results.hubble_parameter(z) for z in z_arr])
    H_arr_invMpc = H_arr / c  # c = 299792.458 km/s

    growth_interp = interp1d(z_arr, D_arr, kind='cubic', fill_value='extrapolate')
    growthf_interp = interp1d(z_arr, f_arr, kind='cubic', fill_value='extrapolate')
    hubble_interp = interp1d(z_arr, H_arr_invMpc, kind='cubic', fill_value='extrapolate')

    return growth_interp, growthf_interp, hubble_interp


# ============================================================
# fixed parameters
# ============================================================
rho_crit = 2.77463e11
AST = 0.3222
aST = 0.707
pST = 0.3
deltaST = 1.686
kappaUV = 1.15e-28
invMpctoinvYear = 3.066067275784999e-07
c = 299792.458
h_HST = 0.7
Omega_m_HST = 0.3
deltaz = 0.5
sigma8_P18 = 0.8102
omega_cdm_P18 = 0.11933
omega_b_P18 = 0.02242
n_s_P18 = 0.9665


# Gaussian quadrature
points, weights = np.polynomial.legendre.leggauss(25)

# Import HST UV LF data
minError = 0.2
DATA_DIR = Path(__file__).resolve().parent
UVLF_PBH = np.loadtxt(DATA_DIR / "UVLF_PBH.txt")
UVLF_PBH[:, 4] = np.array(
    list(map(max, zip(minError * UVLF_PBH[:, 3], UVLF_PBH[:, 4])))
)
UVLF_PBH[:, 5] = np.array(
    list(map(max, zip(minError * UVLF_PBH[:, 3], UVLF_PBH[:, 5])))
)
zs = np.unique(UVLF_PBH[:, 0])

# Beta functions for dust correction
betadata = np.loadtxt(DATA_DIR / "Beta_parameters.txt", unpack=True)
betainterp = PchipInterpolator(betadata[0], betadata[1])
dbetadMUVinterp = PchipInterpolator(betadata[0], betadata[2])

# Undusting the data
dust_corr = []
bin_corr = []
LF_corr = []
for item in UVLF_PBH:
    z, MUV, bin_width = item[0], item[1], item[2]
    new_bin_width = (
        bin_width - AUV(z, MUV + bin_width / 2) + AUV(z, MUV - bin_width / 2)
    )
    dust_corr.append(AUV(z, MUV))
    bin_corr.append(new_bin_width)
    LF_corr.append(bin_width / new_bin_width)

UVLF_PBH[:, 1] -= np.array(dust_corr)
UVLF_PBH[:, 2] = np.array(bin_corr)
UVLF_PBH[:, 3] *= np.array(LF_corr)
UVLF_PBH[:, 4] *= np.array(LF_corr)
UVLF_PBH[:, 5] *= np.array(LF_corr)

# Halo masses
Mhalos = np.geomspace(1e8, 1e14, 1000)

# 全局变量
sigma_cb_interp_func = None
cosmo_initialized = False
growth_interp_global = None
growthf_interp_global = None
hubble_interp_global = None
h_global = None
Omega_m_global = None
Omega_b_global = None
UVLF_data_fixed_global = None
camb_results_global = None


def initialize_fixed_cosmology():
    """
    固定 Planck 2018 宇宙学参数时，只需要初始化一次 CAMB。
    """
    global sigma_cb_interp_func
    global cosmo_initialized
    global growth_interp_global, growthf_interp_global, hubble_interp_global
    global h_global, Omega_m_global, Omega_b_global
    global UVLF_data_fixed_global
    global camb_results_global
    if cosmo_initialized:
        return

    sigma8 = sigma8_P18
    omega_cdm = omega_cdm_P18
    omega_b = omega_b_P18
    n_s = n_s_P18

    pars = camb.CAMBparams()

    pars.set_cosmology(
        H0=None,
        cosmomc_theta=1.041638 / 100.0,
        ombh2=omega_b,
        omch2=omega_cdm,
        tau=0.0544,
        mnu=0.06,
        nnu=3.044,
        num_massive_neutrinos=1,
        theta_H0_range=[10.0, 150.0],
    )

    pars.InitPower.set_params(
        ns=n_s,
        As=2e-9,
        pivot_scalar=0.05,
    )

    z_for_sigma = np.linspace(0, 13, 60)

    pars.set_matter_power(
        redshifts=sorted(z_for_sigma, reverse=True),
        kmax=5.0,
        nonlinear=False,
    )
    pars.WantTransfer = True

    # 第一次 CAMB：校准 sigma8
    results = camb.get_results(pars)
    # 保存固定宇宙学的 CAMB results，供 M_BH(z) 使用
    camb_results_global = results
    h_temp = results.Params.H0 / 100.0
    sigma8_all = results.get_sigmaR(
        8.0 / h_temp,
        hubble_units=False,
        var1='delta_tot',
        var2='delta_tot',
    )
    sigma8_current = float(sigma8_all[-1])

    if sigma8_current <= 0:
        raise RuntimeError("sigma8_current <= 0 in Planck cosmology initialization.")

    As_new = pars.InitPower.As * (sigma8 / sigma8_current) ** 2
    pars.InitPower.set_params(
        ns=n_s,
        As=As_new,
        pivot_scalar=0.05,
    )

    # 第二次 CAMB：最终结果
    results = camb.get_results(pars)

    h_global = results.Params.H0 / 100.0
    Omega_m_global = (
        results.get_Omega('cdm')
        + results.get_Omega('baryon')
        + results.get_Omega('nu')
    )
    Omega_b_global = results.get_Omega('baryon')

    sigma_cb_interp_func = build_sigma_cb_interp(results, z_for_sigma)
    growth_interp_global, growthf_interp_global, hubble_interp_global = build_growth_interps(
        results,
        pars,
    )
    
    UVLF_data_fixed_global = UVLF_PBH.copy()

    for z in zs:
        Vratio = (
            np.power(rcomoving(z + deltaz, Omega_m_HST, h_HST), 3)
            - np.power(rcomoving(z - deltaz, Omega_m_HST, h_HST), 3)
        ) / (
            np.power(rcomoving(z + deltaz, Omega_m_global, h_global), 3)
            - np.power(rcomoving(z - deltaz, Omega_m_global, h_global), 3)
        )

        mask = UVLF_data_fixed_global[:, 0] == z

        UVLF_data_fixed_global[mask, 3:] *= Vratio
        UVLF_data_fixed_global[mask, 1] = UVLF_data_fixed_global[mask, 1] - 5.0 * np.log10(
            rcomoving(z, Omega_m_global, h_global)
            / rcomoving(z, Omega_m_HST, h_HST)
        )
    cosmo_initialized = True
def log_likelihood(param):
    try:
        return log_likelihood_impl(param)
    
    #  捕获所有已知的 CAMB 异常
    except camb.baseconfig.CAMBParamRangeError:
        return -np.inf
    except camb.baseconfig.CAMBError:
        return -np.inf
    except camb.baseconfig.CAMBValueError:
        return -np.inf
    #捕获数值错误
    except (ValueError, FloatingPointError, ZeroDivisionError,
            OverflowError, np.linalg.LinAlgError):
        return -np.inf
    #捕获所有其他异常（最后保底）
    except Exception as e:
        # 记录日志但不崩溃
        import traceback
        with open('camb_errors.log', 'a') as f:
            f.write(f"Unexpected error with params {param}:\n")
            f.write(traceback.format_exc())
            f.write("\n" + "="*60 + "\n")
        return -np.inf


def log_likelihood_impl(param):
    global sigma_cb_interp_func
    (
        logM_PBH,
        logf_PBH,
        loglambda_E68,
        loglambda_E9,
        loglambda_E10,
        loglambda_E11,
        loglambda_E12,
        loglambda_E13p5,
        alphastar,
        betastar,
        epsilonstar_slope,
        epsilonstar_icept,
        Mc_slope,
        Mc_icept,
        sigma_MUV,
        Q,
    ) = param

    initialize_fixed_cosmology()

    h = h_global
    Omega_m = Omega_m_global
    Omega_b = Omega_b_global
    growth_interp = growth_interp_global
    growthf_interp = growthf_interp_global
    hubble_interp = hubble_interp_global

    UVLF_data = UVLF_data_fixed_global
    try:
        MBH_growth_factor_of_z, z_bh_grid, growth_factor_grid = build_MBH_growth_factor_from_camb(
            camb_results=camb_results_global,
            loglambda_E68=loglambda_E68,
            loglambda_E9=loglambda_E9,
            loglambda_E10=loglambda_E10,
            loglambda_E11=loglambda_E11,
            loglambda_E12=loglambda_E12,
            loglambda_E13p5=loglambda_E13p5,
            f_duty=1.0,
            t_Sal=0.045,
            z_merge=20.0,
            z_final=6.0,
            n_points=5000,
        )
    except Exception:
        return -np.inf

    # ============================================================
    # Compute model UV LF
    # ============================================================
    chisq = 0
    for z in zs:
        HMFs = np.array(
            [
                HMF(z, mass, logM_PBH, logf_PBH, h, Omega_b, Omega_m)
                for mass in Mhalos
            ]
        )
        MUV_avs = np.array(
            [
                MUV_from_Mh(
                    growth_interp, growthf_interp, hubble_interp,
                    z, mass, logM_PBH, logf_PBH, h, Omega_b, Omega_m,
                    loglambda_E68, loglambda_E9, loglambda_E10, loglambda_E11, loglambda_E12,
                    alphastar, betastar, epsilonstar_slope, epsilonstar_icept,
                    Mc_slope, Mc_icept, Q, MBH_growth_factor_of_z,
                )
                for mass in Mhalos
            ]
        )

        for item in UVLF_data[UVLF_data[:, 0] == z, 1:]:
            MUV, width, UVLF_ob, UVLF_higherror, UVLF_lowerror = item
            second_integrand = (
                HMFs
                * first_integrand(MUV, width, MUV_avs, sigma_MUV)
                / width
            )
            model_val = simpson(second_integrand, x=Mhalos)
            diff = model_val - UVLF_ob
            if diff > 0:
                chisq += (diff / UVLF_higherror) ** 2
            else:
                chisq += (diff / UVLF_lowerror) ** 2

    return -0.5 * chisq
def uniform_prior_rvs(nwalkers, param_bounds):
    ndim = len(param_bounds)
    pos = np.empty((nwalkers, ndim))
    for i in range(ndim):
        low, high = param_bounds[i]
        pos[:, i] = np.random.uniform(low, high, size=nwalkers)
    return pos

def log_prior(param):
    (
        logM_PBH,
        logf_PBH,
        loglambda_E68,
        loglambda_E9,
        loglambda_E10,
        loglambda_E11,
        loglambda_E12,
        loglambda_E13p5,
        alphastar,
        betastar,
        epsilonstar_slope,
        epsilonstar_icept,
        Mc_slope,
        Mc_icept,
        sigma_MUV,
        Q,
    ) = param
    
    # 1. 基本参数范围检查
    if not (0 < logM_PBH < 14):
        return -np.inf
    if not (-12 < logf_PBH <= 0):
        return -np.inf
    if not (-15 < loglambda_E68 < 2):
        return -np.inf
    if not (-15 < loglambda_E9 < 2):
        return -np.inf
    if not (-15 < loglambda_E10 < 2):
        return -np.inf
    if not (-15 < loglambda_E11 < 2):
        return -np.inf
    if not (-15 < loglambda_E12 < 2):
        return -np.inf
    # 新增的 z >= 13.5 吸积参数
    if not (-15 < loglambda_E13p5 < 2):
        return -np.inf
    if not (-3 < alphastar < 0):
        return -np.inf
    if not (0 < betastar < 3):
        return -np.inf
    if not (-3 < epsilonstar_slope < 3):
        return -np.inf
    if not (-3 < epsilonstar_icept < 3):
        return -np.inf
    if not (-3 < Mc_slope < 3):
        return -np.inf
    if not (7 < Mc_icept < 15):
        return -np.inf
    if not (0.001 < sigma_MUV < 3):
        return -np.inf
    if not (1.5 < Q < 2.5):
        return -np.inf
    
    # 3. 所有检查通过
    return 0
def log_posterior(param):
    lp = log_prior(param)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(param)
# 定义每个参数的均匀分布范围（示例，需按实际情况填写）
param_bounds = [
    (1, 8),           # logM_PBH
    (-12, -1),         # -12 < logf_PBH <= 0
    (-15, 0 ),          # -8 < loglambda_E69 < 2
    (-15, 0 ),
    (-15, 1 ),          # -8 < loglambda_E10 < 2
    (-15, 1 ),          # -8 < loglambda_E11 < 2
    (-15, 1 ),          # -8 < loglambda_E12 < 2
    (-15, 1),           # loglambda_E13p5: 13.5 <= z <= 1000
    (-0.6, -0.2 ),          # -3 < alphastar < 0
    (0.2 , 2.8 ),          # 0 < betastar < 3
    (-2.8, -1 ),          # -3 < epsilonstar_slope < 1
    (-2.8, -1 ),          # -3 < epsilonstar_icept < 0
    (-1.5, 2.5 ),          # -3 < Mc_slope < 3
    (11.5 , 13.5),          # 9< Mc_icept < 13
    (0.001, 1.2),        # 0.001 < sigma_MUV < 3
    (1.6, 2.4),        # 1.5 < Q < 2.5
]
def run_mcmc():
    # 配置 MCMC 参数
    with MPIPool() as pool:
        if not pool.is_master():
            pool.wait()
            return
        ndim = 16     # 参数维度
        nwalkers = 112 # 建议取2-10倍参数维度
        nsteps = 60000   # 链长
        pos = uniform_prior_rvs(nwalkers, param_bounds)
        # 初始化保存文件
        new_backend = emcee.backends.HDFBackend(str(DATA_DIR / 'emcee_UVLF_P18_a16_merged_z20.h5'))
        # 创建进程池并运行采样
        new_sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_posterior, pool=pool, backend=new_backend
        )
        # A fresh backend needs initial positions; an existing backend resumes.
        initial_state = None if new_backend.iteration else pos
        new_sampler.run_mcmc(initial_state, nsteps, progress=True)


if __name__ == "__main__":
    run_mcmc()
