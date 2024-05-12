### SSB.py - a module that computes the next SSB yields given 1, 2, 5, 10-year average bond yields
###
### This module is an implementation of the algorithm published by the Monetary Authority of Singapore (MAS), found at the following URL:
### https://www.mas.gov.sg/-/media/mas/sgs/sgs-announcements-pdf/ssb-pdf/ssb-technical-specifications/ssb-technical-specifications.pdf
### 
### Disclaimer: This author is not affiliated with the Monetary Authority of Singapore. This module is not, shall not be taken as, the
### official implementation of the algorithm used to calculate SSB yields.
### 
### Author: Esther (https://github.com/s-tr)

from __future__ import division

import math
import cvxpy as cp
import numpy as np
from scipy.interpolate import CubicHermiteSpline

# Return an array Y[0..10] denoting the interpolated for each year.
# Slopes determined from the formulas in the "Finite difference" section
# of https://en.wikipedia.org/wiki/Cubic_Hermite_spline
def interpolate(raw_yields):
    Xs = [1,2,5,10]
    ys = raw_yields

    # Determine the slope of the line between adjacent pairs of points
    m = [0,0,0]
    for i in range(3):
        m[i] = (ys[i+1]-ys[i])/(Xs[i+1]-Xs[i])

    # At each point, take the average of the slope of the lines surrounding
    # it. At boundary points, take the slope of the one adjacent line.
    slopes = [0,0,0,0]
    slopes[0] = m[0]
    slopes[1] = (m[0]+m[1])/2
    slopes[2] = (m[1]+m[2])/2
    slopes[3] = m[2]

    spline = CubicHermiteSpline(Xs, ys, slopes)

    # Read off the values of years [1,2,..10]. Value for year 0 can be
    # ignored (it is just there to make the indices nice).
    X_interp = range(11)
    y_interp = spline(X_interp)
    y_interp[0] = 0

    return y_interp

# Return an array DF[0..10] containing the discount factors for years
# 1 to 10. (DF[0] is 1 but is not used)
def generate_discount_factors(Y):
    DF = [1]
    df_sum = 0

    for i in range(1,11):
        DF_i = (1 - Y[i] * df_sum) / (1 + Y[i])
        DF.append(DF_i)
        df_sum += DF_i

    return DF

# Return an array C[0..10] containing the coupon rate of a SSB with
# the zero-arbitrage rule. This set of yields may not satisfy the
# step-up criterion.
def generate_initial_yields(DF):
    sumDFC = 0
    C = [0]

    for i in range(1,11):
        Ci = (1 - sumDFC)/DF[i] - 1
        C.append(Ci)
        sumDFC += DF[i]*C[i]

    return C

# Given a bond with coupon C[0..9] for years 1 to 10, return an array
# C_a[0..9] containing the expected yield if it is redeemed after
# 1..10 years.
def average_yields(C):
    def get_roots(k):
        D = [-1] + C[:k]
        D[-1] += 1
        D = D[::-1]
        return list(np.roots(D))

    def root_filter(root):
        r = np.real(root)
        i = np.imag(root)

        if np.abs(i) <= 1e-6 and 0 <= r and r <= 1:
            return True
        else:
            return False

    roots = []

    for i in range(1,len(C)+1):
        roots_i = get_roots(i)
        root_i = [np.real(r) for r in roots_i if root_filter(r)][0]
        root_i = 1/root_i - 1
        roots.append(root_i)

    return roots


# Given an array containing the 1, 2, 5 and 10 year average bond yields,
# return the predicted Singapore Savings Bond yields for 1..10 years.
# 
# By default, expects and returns values in absolute terms, not percentages.
# e.g. for a value of 3%, input and output will be 0.03.
#
# Pass in `debug=True` to have this function print information on
# intermediate results.
def ssb_solve(raw_yields, debug=False):
    # Utility function.
    def to_percent_string(val):
        return '{:.4f}%'.format(round(val*100,4))
    def print_percents(title, arr):
        if debug:
            print(title + ": " + ", ".join(map(to_percent_string, arr)))

    # 1. Interpolate yields

    Y = interpolate(raw_yields)
    print_percents("Y", Y[1:])

    # 2. Generate discount factors for i = 1..10 based on formula:
    #        1 = DF_1 * Y_i + DF_2 * Y_i + ... DF_i * Y_i + DF_i * 1
    #     :: DF_i = (1 - Y_i * (DF_1 + DF_2 + ... + DF_{i-1})) / (1 + Y_i)

    DF = generate_discount_factors(Y)
    print_percents("D", DF[1:])

    # 3. Generate initial yields which may not be step-up
    # This step is not strictly necessary, but we include it as it is part of
    # the official spec.
    C0 = generate_initial_yields(DF)
    print_percents("C0", C0[1:])
    print_percents("Y0", average_yields(C0[1:]))

    # 4. Define coupon yields C_i for i = 1..10.

    C = [0]*11
    for i in range(1,11):
        C[i] = cp.Variable(1)

    # 5. Define the value e_i for i = 1..10:
    #        e_i = DF_1 * C_1 + DF_2 * C_2 + ... DF_i * C_i + DF_i * 1 - 1
    # e_i can be understood as the spread between holding an SSB for i years
    # instead of a government bond maturing in i years.
    # e_i > 0 when holding SSB gives more return, and vice versa.

    E = [0] * 11
    for i in range(1,11):
        e_i = -1
        for j in range(1, i+1):
            e_i += DF[j] * C[j]
        e_i += DF[i]
        E[i] = e_i

    # 6. Constraints
    constraints = []

    # a. coupons are non-decreasing
    for i in range(10):
        constraints.append(C[i] <= C[i+1])

    # b. yields must not exceed yield of equivalent government bonds
    for i in range(1,10):
        constraints.append(E[i] <= 0)

    # c. 10-year yield must be equal to 10-year government bond yield
    constraints.append(E[10] == 0)

    # 7. Find the values of C which minimise sum(e_i^2) for i=1..10.

    objective_function = 0
    for i in range(1,11):
        objective_function += E[i]**2
    objective = cp.Minimize(objective_function)

    opt = cp.Problem(objective, constraints).solve(eps_abs=1e-6,eps_rel=1e-6)

    # 8. Read off the values of C.
    C_opt = [C[i].value[0] for i in range(1,11)]
    print_percents("C*", C_opt)
    print_percents("Y*", average_yields(C_opt))

    return C_opt

# Given an array containing the 1, 2, 5 and 10 year average bond yields,
# return the predicted Singapore Savings Bond yields for 1..10 years.
# 
# This function expects and returns values in percentage terms,
# e.g. for a value of 3%, input and output will be 3.
#
# By default, results are rounded to 2 decimal places. To get the raw
# results without rounding, pass `rounding=False`.
#
# Pass in `debug=True` to have this function print information on
# intermediate results.
def ssb(benchmark_yields, debug=False, rounding=True):
    benchmark_yields = [x/100 for x in benchmark_yields]

    raw_result = ssb_solve(benchmark_yields, debug=debug)
    if rounding:
        result_percents = [round(100*x + 0.0000001, 2) for x in raw_result]
        if debug:
            resid = [round(10000*(p-100*x)) for p,x in zip(result_percents,raw_result)]
            print("Rounding: " + str(resid))
    else:
        result_percents = [100*x for x in raw_result]

    return result_percents

