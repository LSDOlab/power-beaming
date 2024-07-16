import numpy as np

def standard_atmosphere(alt):
    R = 287 # [J/kg-K] specific gas constant for air
    if alt > 25000:
        T = -131.21 + 0.00299 * alt
        rho = (2.488 * (((T + 273.1) / 216.6) ** -11.388)) / (0.2869 * (T + 273.1))
    elif 11000 < alt <= 25000:
        T = -56.46
        rho = (22.65 * np.exp(1.73 - 0.000157 * alt)) / (0.2869 * (T + 273.1))
    elif alt <= 11000:
        T = 15.04 - 0.00649 * alt
        rho = (101.29 * (((T + 273.1) / 288.08) ** 5.256)) / (0.2869 * (T + 273.1))

    T = T + 273.15 # result in K 
    p = rho*R*T / 100 # result in hPa    

    return p, T