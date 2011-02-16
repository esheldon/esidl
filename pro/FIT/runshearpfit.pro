PRO runshearpfit

dir='/sdss4/data1/esheldon/CLUSTER/'
f=strarr(3)
f[0] = dir+'zcluster_752_756_g_N1.dat'
f[1] = dir+'zcluster_752_756_r_N1.dat'
f[2] = dir+'zcluster_752_756_i_N1.dat'

; initial guesses
asig = [40., -.7]
ashear = [.005, -.7]
type = 2
names = ['g','r','i']
yrange = [-50, 200]
title='26 Clusters'

shearpfit,type,f,ashear,names=names,title=title
key=get_kbrd(1)

shearpfit,type,f,asig,tag='SIGMA',/model,names=names,yrange=yrange,title=title

return
END 
