function mpfit_schechter_gauss, m, p

    phistar = p[0]
    mstar = p[1]
    alpha = p[2]

    norm = p[3]
    mean= p[4]
    sig = p[5]

    func = $
        0.4*alog(10.)*phistar*10.^(-0.4*(m-mstar)*(alpha+1.))* $
        exp(-10.^(-0.4*(m-mstar))) + norm*exp( (m-mean)^2/2.0/sig^2)
    return, func

end
