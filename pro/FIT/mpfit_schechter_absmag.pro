function mpfit_schechter_absmag, m, p

    phistar = p[0]
    mstar = p[1]
    alpha = p[2]
    func = $
        0.4*alog(10.)*phistar*10.^(-0.4*(m-mstar)*(alpha+1.))* $
        exp(-10.^(-0.4*(m-mstar)))
    return, func

end
