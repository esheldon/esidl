function nfw_ds_plus_lin,r,p

linear_ds,r,dslin

return,nfw_delta_sigma(r,p[0:1])+p[2]*dslin
end
