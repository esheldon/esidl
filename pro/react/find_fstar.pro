PRO find_fstar,fstar, z, sigmasq

m = n_elements(z)
g = (z^2.0 - sigmasq)/z^2.0
w = z^2/double(m)

fstar = antitonic(g,w,0)

RETURN
END
