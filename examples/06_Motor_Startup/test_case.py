"""PYPOWER power flow data for 2 bus, 1 generator case.
"""

from numpy import array

def test_case():
    """PYPOWER power flow data for 2 bus, 1 generator case.
    """
    ppc = {"version": '2'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    ppc["bus"] = array([
        [1, 3, 0,  0, 0, 0, 1, 0.95, 0, 345, 1, 1.1, 0.9],
        [2, 1, 5,  2, 0, 0, 1, 1, 0, 345, 1, 1.1, 0.9]
    ])

    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
    ppc["gen"] = array([
        [1, 0,   0, 300, -300, 1.0, 100, 1, 250, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ])

    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["branch"] = array([
        [1, 2, 0.01,   0.0576, 0, 250, 250, 250, 0, 0, 1, -360, 360],
        [1, 2, 0.01,   0.085,  0, 250, 250, 250, 0, 0, 1, -360, 360]
    ])

    return ppc
