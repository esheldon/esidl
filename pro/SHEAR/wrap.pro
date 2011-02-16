PRO wrap,scat=scat

lcat=mrdfits('/sdss4/data1/esheldon/CORRECTED/run752_756_lensgal_r_overlap.fit',1)


;runkappa_map,lcat,slength=10.,rfac=15.,gridsize=100.,stepfac=2.,scat=scat,random=30,/silent,/sum,/surface,/write,/noprompt

;runkappa_map,lcat,slength=5.,rfac=15.,gridsize=100.,stepfac=2.,scat=scat,random=30,/silent,/sum,/surface,/write,/noprompt

;runkappa_map,lcat,slength=10.,rfac=10.,gridsize=100.,stepfac=2.,scat=scat,random=60,/silent,/sum,/noprompt,/ret

;runkappa_map,lcat,slength=10.,rfac=10.,gridsize=100.,stepfac=2.,scat=scat,/silent,/sum,/surface,/write,/noprompt,/allign,noise=0.605711

;w = where(lcat.r LE .35)

;runkappa_map,lcat[w],slength=10.,rfac=10.,gridsize=100.,stepfac=2.,scat=scat,/silent,/sum,/surface,/write,/noprompt,/allign,noise=0.605711

runkappa_map,lcat,slength=120.,rfac=10.,stepfac=2.,scat=scat,random=150,/silent,/sum,/surface,/noprompt,/ret,neach=1

return
END 
