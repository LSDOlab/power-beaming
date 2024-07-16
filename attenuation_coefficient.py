import numpy as np

# From ITU-R P 676
# This code is valid for 1 to 1000 GHz (300 micrometers - 0.3 meters)

def attenuation_coefficient(f, alt, p_tot, T):

    rho_H2O_0 = 7.5 # [g/m^3] Correct units?

    rho_H2O = rho_H2O_0*np.exp(-alt/2000)

    e = rho_H2O*T/216.7 # water vapor partial pressure. Units in paper are confusing, T in K or C?

    p_air = p_tot - e

    theta = 300/T

    OX_data = "v12_lines_oxygen.txt"

    H2O_data = "v12_lines_water_vapour.txt"

    OX_table = np.loadtxt(OX_data, dtype=float, skiprows=1,delimiter=',')

    H2O_table = np.loadtxt(H2O_data, dtype=float, skiprows=1,delimiter=',')

    f0s_OX = OX_table[:,0]
    a1s = OX_table[:,1]
    a2s = OX_table[:,2]
    a3s = OX_table[:,3]
    a4s = OX_table[:,4]
    a5s = OX_table[:,5]
    a6s = OX_table[:,6]
    Sis_OX = [0] * len(f0s_OX)
    Deltaf_OX = [0] * len(f0s_OX)
    delta_OX = [0] * len(f0s_OX)
    NppD_OX = 0
    Fis_OX = [0] * len(f0s_OX)

    f0s_H2O = H2O_table[:,0]
    b1s = H2O_table[:,1]
    b2s = H2O_table[:,2]
    b3s = H2O_table[:,3]
    b4s = H2O_table[:,4]
    b5s = H2O_table[:,5]
    b6s = H2O_table[:,6]
    Sis_H2O = [0] * len(f0s_H2O)
    Deltaf_H2O = [0] * len(f0s_H2O)
    Fis_H2O = [0] * len(f0s_H2O)

    for i in range(len(f0s_OX)):
        Sis_OX[i] = a1s[i] * 10**-7 * p_air * theta**3 * np.exp(a2s[i] * (1 - theta))
        Deltaf_OX[i] = a3s[i] * 10**-4 * (p_air * theta**(0.8 - a4s[i]) + 1.1 * e * theta)
        Deltaf_OX[i] = np.sqrt(Deltaf_OX[i]**2 + 2.25 * 10**-6)
        delta_OX[i] = (a5s[i] + a6s[i]*theta) * 10**-4 * (p_air + e) * theta**0.8

    if f < 10 or f > 100:
        d = 5.6 * 10**-4 * (p_air + e) * theta**0.8
        NppD_OX = f * p_air * theta**2 * ((6.14 * 10**-5)/(d*(1 + (f/d)**2)) + (1.4 * 10**-12 * p_air * theta**1.5)/(1 + 1.9 * 10**-5 * f**1.5))

    for l in range(len(Fis_OX)):
        Fis_OX[l] = f/f0s_OX[l] * ((Deltaf_OX[l] - delta_OX[l] * (f0s_OX[l] - f))/((f0s_OX[l] - f)**2 + Deltaf_OX[l]**2) + (Deltaf_OX[l] - delta_OX[l] * (f0s_OX[l] + f))/((f0s_OX[l] + f)**2 + Deltaf_OX[l]**2))

    for j in range(len(f0s_H2O)):
        Sis_H2O[j] = b1s[j] * 10**-1 * e * theta**3.5 * np.exp(b2s[j]*(1-theta))
        Deltaf_H2O[j] = b3s[j] * 10**-4 * (p_air * theta**b4s[j] + b5s[j] * e * theta**b6s[j])
        Deltaf_H2O[j] = 0.535*Deltaf_H2O[j] + np.sqrt(0.217 * Deltaf_H2O[j]**2 + 2.136 * 10**-12 * f0s_H2O[j]**2/theta)

    for m in range(len(Fis_H2O)):
        Fis_H2O[m] = f/f0s_H2O[m] * ((Deltaf_H2O[m])/((f0s_H2O[m] - f)**2 + Deltaf_H2O[m]**2) + (Deltaf_H2O[m])/((f0s_H2O[m] + f)**2 + Deltaf_H2O[m]**2))

    Npp_H2O = np.dot(Sis_H2O,np.transpose(Fis_H2O))

    Npp_OX = np.dot(Sis_OX,np.transpose(Fis_OX)) + NppD_OX

    gamma = 0.1820*f*(Npp_OX + Npp_H2O) * 0.001

    return gamma