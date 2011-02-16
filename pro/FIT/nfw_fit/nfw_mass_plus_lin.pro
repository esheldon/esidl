function nfw_mass_plus_lin,r,p
rhoc=0.27752
j3=interpol(*!linear_j3_cache,*!linear_r_cache,r,/spline)

return,nfw_mass(r,p[0:1])+p[2]*rhoc*j3
end
