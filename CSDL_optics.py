import numpy as np
import csdl_alpha as csdl

def CSDL_optics(d_t, d_r, lam, n, dist, theta, alt_t, alt_r, alt_steps):

    '''This code takes transmitter diameter, receiver diameter, beam wavelength,
    refractive index, and beam divergence angle as inputs, and calculates
    whether or not the beam diameter exceeds the receiver diameter at the
    given distance. If so, it calculates the resultant energy losses.'''

    # Calculate beam waist and Rayleigh range
    w0 = lam / (np.pi * n * theta)  # [m] (beam waist)
    z_R = np.pi * w0**2 * n / lam # (distance at which intensity is half of its peak at the waist)

    z_R_max = np.pi*(d_r/2)**2*n/lam # where beam waist equals receiver width. 
    # Any higher than this and you could NEVER receive 100% of the energy cone

    # if z_R > z_R_max:
    #     print('Warning: beam waist exceeds receiver width. 100% energy capture cannot be achieved.')

    # Solve for the distance at which the beam diameter equals the transmitter diameter
    dist_t = csdl.ImplicitVariable(name='dist_t', value=dist/2)

    # if 2*w0 > d_t:
    #     print('Configuration not possible; beam waist larger than transmitter diameter. Adjust parameters and rerun.')
    #     raise SystemExit
    
    ''' dist_t is the distance from the waist (in either direction) at which the beam width is
    # equal to the transmitter width. Practically this gives the distance
    # between the waist and the transmitter '''
    
    residual_1 = d_t / 2 - w0 * csdl.sqrt(1 + (dist_t/ z_R)**2)
    residual_1.add_name('residual_1')

    solver = csdl.nonlinear_solvers.Newton('nlsolver_z')
    solver.add_state(dist_t, residual_1) # Specify that the residual of the state x1 is residual_1 
    
    solver.run()

    # Need to make dist_t positive if for some reason the solver gives its negative conjugate.
    dist_t = (dist_t**2)**0.5

    # print('dist_t: ',dist_t.value)
    # print('residual: ',residual_1.value)

    # Calculate distance from waist to receiver
    dist_r = ((dist-dist_t)**2)**0.5

    ''' beam width at receiver. This value should be less than the width 
    of the reciever so that the full cone of energy is captured '''
    w_r = w0 * csdl.sqrt(1 + (dist_r / z_R)**2)

    eta_ts = [0] * len(alt_steps)
    j = 0

    for i in alt_steps:
        eta_ts[j] =  d_r/2 / (w0 * csdl.sqrt(1 + ( ((i - dist_t)**2)**0.5 / z_R)**2))
        j = j + 1

    wtf = 0
    for i in eta_ts:
        eta_ts[wtf] = eta_ts[wtf] - 1/(1+csdl.exp(-1000*(eta_ts[wtf]-1))) * (eta_ts[wtf]-1)
        # print(eta_ts[wtf].value)
        wtf = wtf + 1

    loss = (1 - eta_ts[-1])*100 # fractional amount of power wasted from 0-1

    # Print the result
    # if eta_t < 1:
    #     print('The beam diameter is larger than the receiver diameter at the given distance.\n{loss}% of the available energy is lost.')
    #     print('Consider adjusting beam divergence or increasing the diameter of the transmitter or receiver.')
    # else:
    #     print('The receiver is able to capture 100% of the available energy.')
    
    print('The receiver is able to capture',  eta_ts[-1].value * 100,'% of the available energy.')

    return eta_ts


# Giving up on this but in the future I'd like to be able to actively see loss due to optical diffraction
    # dist_steps = np.linspace(alt_t,dist,int((alt_r-alt_t)/stepsize))

    # dists_from_waist = [0] * len(dist_steps)
    # ws = [0] * len(dist_steps)
    # for i in range(len(dist_steps)):
    #     dists_from_waist[i] = dist_t - dist_steps[i]
    #     ws = w0 * np.sqrt(1-(dists_from_waist[i]/z_R)**2)

    # if dist_t > dist and 