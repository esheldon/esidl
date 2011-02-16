pro define_nfw_200,z=z,Omega_m=Omega_m,Omega_l=Omega_l
;define params for NFW M200 stuff
;assumes simple flat cosmology
;assumes z=0.25 by default, close to cluster mean

if n_elements(z) eq 0 then z=0.25
if n_elements(Omega_m) eq 0 then Omega_m=0.27

defsysv,'!nfw_200',{z:z,Omega_m:Omega_m}
return
end
