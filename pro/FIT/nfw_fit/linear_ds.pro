pro linear_ds,r,ds
;computes the comoving signal and then
;redshifts it properly

z=!nfw_200.z
Omega_m=!nfw_200.Omega_m

linear_ds_com,r*(1+z),ds
;growth_factor,1.0/(1+z),D,Omega_m
;leave out the growth factor
D=1.0
ds=ds*D^2*(1.+z)^2

;now the data will constrain B=Omega_m*sigma_8^2*D(z)^2*b(M,z)

return
end


