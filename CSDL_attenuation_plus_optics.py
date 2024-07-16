# to do:
# replace fixed refractive index with ITU-R rec

import csdl_alpha as csdl
import numpy as np
import matplotlib.pyplot as plt
from attenuation_coefficient import *
from standard_atmosphere import *
from CSDL_optics import *

# variables represent inputs to the model with fixed or changing values.
lam = 0.001  # [m] (beam wavelength)
I_t = 10000 # input power
alt_t = 0 # [m] altitude of transmitter
alt_r = 5000 # [m] altitude of receiver
n = 1.0003  # (refractive index) REPLACE THIS WITH ITU-R REC
stepsize = 100 # [m] determines how often absorption is re-calculated. Lower number means more accurate results but more computation

c = 299792458 # [m/s]
f =  c/lam * 10**-9 # [GHz]

# You should need to specify only one of these
slant_angle = 90 # slant angle, measured from line tangent to ground
# horizontal_distance = 1000 

if 'slant_angle' in locals() and 'horizontal_distance' in locals():
    print('Having both a slant angle and a horizontal distance is redundant. Comment out one and leave the other.')
    raise SystemExit

if f > 1000 or f < 1:
    print('''The current version of this code performs its atmospheric attenuation calculations based on ITU-R P.676-13,
        which is only valid for electromagnetic frequencies ranging from 1 to 1000 GHz (30 cm to 3000 micrometers).
        Your chosen frequency is outside of this range, so attenuation estimates will not be accurate.)''')

# loss due to attenuation (absorption and scattering)
if alt_t == alt_r:
    if 'horizontal_distance' in locals():
        steps = np.ceil(horizontal_distance/stepsize) + 1 # may need to revisit this
        alt_steps = [alt_t] * steps.astype(np.int64)
    else:
        print('Since this is a horizontal path, horizontal distance must be specified. See "Inputs" at top of code.')
        raise SystemExit
else:
    alt_steps = list(range(alt_t, alt_r, stepsize))
    alt_steps.append(alt_r)
# initialize arrays
ps = [0]*len(alt_steps)
Ts = [0]*len(alt_steps)
gammas = [0]*len(alt_steps)
i = 0

for x in alt_steps:
    ps[i], Ts[i] = standard_atmosphere(x)
    gammas[i] = attenuation_coefficient(f, alt_steps[i], ps[i], Ts[i]) # result in [dB/m]
    # gammas[i] = 10**(-gammas[i]/10)/1000 # result in [1/m]
    i = i + 1

Is = [0]*len(alt_steps)
Is[0] = I_t
j = 1

for g in gammas:
    Is[j] = Is[j-1]*10**(-gammas[j]/10 * stepsize)
    j = j + 1
    if j == len(alt_steps):
        break

# loss due to optical diffraction
if 'slant_angle' in locals():
    dist = (alt_r-alt_t)/np.sin(slant_angle*np.pi/180)
else:
    dist = horizontal_distance

recorder = csdl.Recorder(inline=True)
recorder.start()

# Values that will be optimized by CSDL
d_t = csdl.Variable(value=5,name='Transmitter Diameter')  # [m] (transmitter diameter) WILL BE A CSDL VARIABLE
d_r = csdl.Variable(value=2,name='Receiver Diameter')  # [m] (receiver diameter) WILL BE A CSDL VARIABLE
theta = csdl.Variable(value=4e-4,name='Beam Divergence')  # [rad] (beam divergence, dependent on the lens chosen)

d_t.set_as_design_variable(lower=0.1,upper=10)
d_r.set_as_design_variable(lower=0.1,upper=10)
theta.set_as_design_variable(lower=1e-5,upper=1e-3)

eta_ts = CSDL_optics(d_t, d_r, lam, n, dist, theta, alt_t, alt_r, alt_steps)

# THIS WILL NOT WORK IN OPTIMIZATION STEP
# for i in eta_ts:
#     if i.value > 1:
#         i.value = 1
#     print(i.value)

I_rs = [0] * len(alt_steps)

# for x in range(len(eta_ts)):
#     print(eta_ts[x].value)

j = 0
for i in alt_steps:
    print(eta_ts[j].value)
    I_rs[j] = Is[j] * eta_ts[j]
    # print(I_rs[j])
    j = j + 1

# print(eta_ts[j-1].value)
# print(Is[j-1])
# print("If " + str(I_t) + " Watts of energy is beamed, " + I_rs[-1].value + " W will be received.")

recorder.stop()

# Is_list = [0] * len(Is)
I_rs_list = [0] * len(Is)

for i in range(len(Is)):
    # Is_list[i] = Is[i].value
    I_rs_list[i] = I_rs[i].value

if alt_t == alt_r:
    plt.plot(list(range(0, horizontal_distance + stepsize, stepsize)),Is)
    plt.xlabel("Horizontal Distance [m]") 
else:
    plt.plot(alt_steps,Is,'-', label='Received power not considering optics')
    plt.plot(alt_steps,I_rs_list,'--',label='Received power considering optics')
    plt.xlabel("Altitude [m]")

plt.ylabel("Received Power")
plt.legend()
plt.show()

# recorder.active_graph.visualize()