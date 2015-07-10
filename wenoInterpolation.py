def JSWENO_interpolation(sp0,sp1,sp2,sp3,sp4):
    '''
    Input:
    ------
            Five stencil points
    Output:
    -------
            Interpolated value at x+.5
    '''

    epsilon = 10**-12  # to avoid dividing by zero

    # polynomial approximation of X for three stencils
    p0 = (1./3.)*sp0 - (7./6.)*sp1 + (11./6.)*sp2
    p1 = -(1./6.)*sp1 + (5./6.)*sp2 + (1./3.)*sp3
    p2 = (1./3.)*sp2 + (5./6.)*sp3 - (1./6.)*sp4

    gamma = [.1, .6, .3]  # ideal (linear) weights

    # nonlinear weights (omegas):
    # The weight is determined by the smoothness of each stencil

    # betas are smoothness indicators for each stencil
    # they come from a complicated formula (cf Borges et al. 2008)
    # basically they measure the regularity (or variation) of u over each stencil
    beta0 = (13./12.)*(sp0 - 2*sp1 + sp2)**2 + .25*(sp0 - 4*sp1 + 3*sp2)**2
    beta1 = (13./12.)*(sp1 - 2*sp2 + sp3)**2 + .25*(sp1 - sp3)**2
    beta2 = (13./12.)*(sp2 - 2*sp3 + sp4)**2 + .25*(3*sp2 - 4*sp3 + sp4)**2


    # alphas are the unnormalized weights
    alpha0 = gamma[0] / (epsilon + beta0)**2
    alpha1 = gamma[1] / (epsilon + beta1)**2
    alpha2 = gamma[2] / (epsilon + beta2)**2

    alpha_sum = alpha0 + alpha1 + alpha2

    # normalized nonlinear weights of Jiang and Shu
    w0JS = alpha0 / alpha_sum
    w1JS = alpha1 / alpha_sum
    w2JS = alpha2 / alpha_sum

    u_half = w0JS * p0 + w1JS * p1 + w2JS * p2

    return u_half

def WENOM_interpolation(sp0,sp1,sp2,sp3,sp4):
    '''
    Input:
    ------
            Five stencil points
    Output:
    -------
            Interpolated value at x+.5
    '''

    epsilon = 10**-12  # to avoid dividing by zero

    # polynomial approximation of X for three stencils
    p0 = (1./3.)*sp0 - (7./6.)*sp1 + (11./6.)*sp2
    p1 = -(1./6.)*sp1 + (5./6.)*sp2 + (1./3.)*sp3
    p2 = (1./3.)*sp2 + (5./6.)*sp3 - (1./6.)*sp4

    gamma = [.1, .6, .3]  # ideal (linear) weights

    # nonlinear weights (omegas):
    # The weight is determined by the smoothness of each stencil

    # betas are smoothness indicators for each stencil
    # they come from a complicated formula (cf Borges et al. 2008)
    # basically they measure the regularity (or variation) of u over each stencil
    beta0 = (13./12.)*(sp0 - 2*sp1 + sp2)**2 + .25*(sp0 - 4*sp1 + 3*sp2)**2
    beta1 = (13./12.)*(sp1 - 2*sp2 + sp3)**2 + .25*(sp1 - sp3)**2
    beta2 = (13./12.)*(sp2 - 2*sp3 + sp4)**2 + .25*(3*sp2 - 4*sp3 + sp4)**2


    # alphas are the unnormalized weights
    alpha0 = gamma[0] / (epsilon + beta0)**2
    alpha1 = gamma[1] / (epsilon + beta1)**2
    alpha2 = gamma[2] / (epsilon + beta2)**2

    alpha_sum = alpha0 + alpha1 + alpha2

    # normalized nonlinear weights of Jiang and Shu
    w0JS = alpha0 / alpha_sum
    w1JS = alpha1 / alpha_sum
    w2JS = alpha2 / alpha_sum

    # Henrick Mapping:
    def g(k, w):
        return (w * (gamma[k] + gamma[k]**2 - 3*gamma[k]*w + w**2))/(gamma[k]**2 + w*(1-2*gamma[k]))

    alphas_0 = g(0, w0JS)
    alphas_1 = g(1, w1JS)
    alphas_2 = g(2, w2JS)

    alphas = [alphas_0, alphas_1, alphas_2]

    # mapped weights
    w0m = alphas_0/sum(alphas)
    w1m = alphas_1/sum(alphas)
    w2m = alphas_2/sum(alphas)

    u_half = w0m * p0 + w1m * p1 + w2m * p2

    return u_half

def WENOZ_interpolation(sp0,sp1,sp2,sp3,sp4):
    '''
    New smoothness indicators based on previous JS betas
    The new smoothness indicators are higher order
    Tau = | Beta0 - Beta2 |
    Betaz_k = (Betak + epsilon) / (Betak + Tauk + epsilon)
    '''

    epsilon = 10**-12  # to avoid dividing by zero

    # polynomial approximation of X for three stencils
    p0 = (1./3.)*sp0 - (7./6.)*sp1 + (11./6.)*sp2
    p1 = -(1./6.)*sp1 + (5./6.)*sp2 + (1./3.)*sp3
    p2 = (1./3.)*sp2 + (5./6.)*sp3 - (1./6.)*sp4

    gamma = [.1, .6, .3]  # ideal (linear) weights

    # nonlinear weights (omegas):
    # The weight is determined by the smoothness of each stencil

    # betas are smoothness indicators for each stencil
    # they come from a complicated formula (cf Borges et al. 2008)
    # basically they measure the regularity (or variation) of u over each stencil
    beta0 = (13./12.)*(sp0 - 2*sp1 + sp2)**2 + .25*(sp0 - 4*sp1 + 3*sp2)**2
    beta1 = (13./12.)*(sp1 - 2*sp2 + sp3)**2 + .25*(sp1 - sp3)**2
    beta2 = (13./12.)*(sp2 - 2*sp3 + sp4)**2 + .25*(3*sp2 - 4*sp3 + sp4)**2

    # Borges WENOZ smoothness indicators:
    tau = abs(beta2 - beta0)

    alpha0z = gamma[0] * ( 1 + (tau / (beta0 + epsilon))**2)
    alpha1z = gamma[1] * ( 1 + (tau / (beta1 + epsilon))**2)
    alpha2z = gamma[2] * ( 1 + (tau / (beta2 + epsilon))**2)

    alphaz_sum = alpha0z + alpha1z + alpha2z

    # normalized nonlinear weights
    w0z = alpha0z / alphaz_sum
    w1z = alpha1z / alphaz_sum
    w2z = alpha2z / alphaz_sum


    u_half = w0z * p0 + w1z * p1 + w2z * p2

    return u_half

def WENOSG_interpolation(sp0,sp1,sp2,sp3,sp4):
    """
    Doesnt work for whatever reason
    """

    epsilon = 10**-12

    p0 = (1./3.)*sp0 - (7./6.)*sp1 + (11./6.)*sp2
    p1 = -(1./6.)*sp1 + (5./6.)*sp2 + (1./3.)*sp3
    p2 = (1./3.)*sp2 + (5./6.)*sp3 - (1./6.)*sp4

    gamma = [.1, .6, .3]

    beta0 = (13./12.)*(sp0 - 2*sp1 + sp2)**2 + .25*(sp0 - 4*sp1 + 3*sp2)**2
    beta1 = (13./12.)*(sp1 - 2*sp2 + sp3)**2 + .25*(sp1 - sp3)**2
    beta2 = (13./12.)*(sp2 - 2*sp3 + sp4)**2 + .25*(3*sp2 - 4*sp3 + sp4)**2

    c00 = .25
    c10 = .75
    c01 = .5
    c11 = .5

    tau0 = abs(beta1 - beta0)
    tau1 = abs(beta2 - beta1)

    alpha00 = c00*(1 + (tau0 / (beta0 + epsilon))**2)
    alpha10 = c10*(1 + (tau0 / (beta1 + epsilon))**2)
    alpha01 = c01*(1 + (tau1 / (beta1 + epsilon))**2)
    alpha11 = c11*(1 + (tau1 / (beta2 + epsilon))**2)

    alphakl_sum = alpha00 + alpha10 + alpha01 + alpha11

    psi00 = alpha00 / alphakl_sum
    psi10 = alpha11 / alphakl_sum
    psi01 = alpha01 / alphakl_sum
    psi11 = alpha11 / alphakl_sum

    def g(k, w):
        return (w * (gamma[k] + gamma[k]**2 - 3*gamma[k]*w + w**2))/(gamma[k]**2 + w*(1-2*gamma[k]))

    alphastar00 = g(0, psi00)
    alphastar10 = g(1, psi10)
    alphastar01 = g(0, psi01)
    alphastar11 = g(1, psi11)

    alphastar_sum = alphastar00 + alphastar10 + alphastar01 + alphastar11

    w00 = alphastar00 / alphastar_sum
    w10 = alphastar10 / alphastar_sum
    w01 = alphastar01 / alphastar_sum
    w11 = alphastar11 / alphastar_sum

    h0 = w00*p0 + w10*p1
    h1 = w01*p1 + w11*p2

    tau = abs(beta2 - beta0)

    c0 = .4
    c1 = .6

    alpha0 = c0 * (1 + (tau / (beta0 + epsilon))**2)
    alpha1 = c1 * (1 + (tau / (beta2 + epsilon))**2)

    alpha_sum = alpha0 + alpha1

    psi0 = alpha0 / alpha_sum
    psi1 = alpha1 / alpha_sum

    psi_sum = g(0, psi0) + g(1, psi1)

    w0SG = g(0, psi0) / psi_sum
    w1SG = g(1, psi1) / psi_sum

    u_half = w0SG * h0 + w1SG * h1

    return u_half
