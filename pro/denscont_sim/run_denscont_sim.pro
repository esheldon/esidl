PRO run_denscont_sim

  central_sigmav = 120d
  fac = 10d
  neighbor_sigmav = 10d*central_sigmav
  neighbor_pos = -500d          ;kpc
  Npts = 2500L

  base = 'denscont_sim_cent120_neigh10times_500kpc'
  psfile = '/home/esheldon/idl.lib/denscont_sim/'+base+'.ps'
  savefile = '/home/esheldon/idl.lib/denscont_sim/'+base+'.sav'

  begplot, name = psfile, /color
  denscont_sim, central_sigmav, neighbor_sigmav, neighbor_pos, Npts, $
                rmax, $
                density_area_central, density_circ_central, density_contrast_central,$
                density_area_neigh, density_circ_neigh, density_contrast_neigh,$
                density_area, density_circ, density_contrast

  density_cont_sim_plot, rmax, $
                         density_area_central, density_circ_central, density_contrast_central,$
                         density_area_neigh, density_circ_neigh, density_contrast_neigh,$
                         density_area, density_circ, density_contrast, /loglog
  endplot

  save, rmax, $
        density_area_central, density_circ_central, density_contrast_central,$
        density_area_neigh, density_circ_neigh, density_contrast_neigh,$
        density_area, density_circ, density_contrast, $
        file = savefile
  

END 
