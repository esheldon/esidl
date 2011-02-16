function nfw_sigma_core_sys,r
;return NFW sigma(R)
;with a core
;requires systems variables set
;!nfw_sigma_core_sys_core
;!nfw_sigma_core_sys_rvir
;!nfw_sigma_core_sys_c

core=!nfw_sigma_core_sys_core
x=sqrt(r^2+ core^2)

rvir=!nfw_sigma_core_sys_rvir
c=!nfw_sigma_core_sys_c

p=[rvir,c]
return,nfw_sigma(x,p)
end
