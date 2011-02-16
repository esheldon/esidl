PRO runsimpwtheta

indir='/sdss4/data1/esheldon/CORRECTED/'

lcat=mrdfits(indir+'run752_756_lensgal_g_overlap.fit',1)
clr=1
voronoi_simpwtheta, lcat, clr

lcat=mrdfits(indir+'run752_756_lensgal_r_overlap.fit',1)
clr=2
voronoi_simpwtheta, lcat, clr

lcat=mrdfits(indir+'run752_756_lensgal_i_overlap.fit',1)
clr=3
voronoi_simpwtheta, lcat, clr

return
END 
