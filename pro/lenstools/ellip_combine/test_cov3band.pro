PRO test_cov3band_gete, struct, rcut, errcut, $
               e1corr, e2corr, e1e1err, e1e2err, e2e2err, w

  defval = -9999.

  nobj = n_elements(struct)

  e1corr = replicate(defval, 3, nobj)
  e2corr = e1corr
  e1e1err = e1corr
  e1e2err = e1corr
  e2e2err = e2corr

  corr = struct.m_r[1:3]

  g=0
  r=1
  i=2

  ;; correct r-band e1/e2
  e1corr[r,*]  = struct.m_e1_corr[2]
  e2corr[r,*]  = struct.m_e2_corr[2]
  e1e1err[r,*] = struct.m_e1e1err[2]
  e1e2err[r,*] = struct.m_e1e2err[2]
  e2e2err[r,*] = struct.m_e2e2err[2]

  e1corr[r,*] = e1corr[r,*]/(1. - corr[r,*])
  e2corr[r,*] = e2corr[r,*]/(1. - corr[r,*])
  e1e1err[r,*] = e1e1err[r,*]/(1. - corr[r,*])
  e1e2err[r,*] = e1e2err[r,*]/(1. - corr[r,*])
  e2e2err[r,*] = e2e2err[r,*]/(1. - corr[r,*])

  ;; look for good i-band measurements

  w=where(struct.m_r[3] GT 0.0 AND struct.m_r[3] LT rcut, niband)

  w2=where(struct[w].m_e1e1err[3] LT errcut AND $
           struct[w].m_e2e2err[3] LT errcut, niband)
  w=w[w2]

  ;; correct i-band stuff
  e1corr[i,w]  = struct[w].m_e1_corr[3]
  e2corr[i,w]  = struct[w].m_e2_corr[3]
  e1e1err[i,w] = struct[w].m_e1e1err[3]
  e1e2err[i,w] = struct[w].m_e1e2err[3]
  e2e2err[i,w] = struct[w].m_e2e2err[3]

  e1corr[i,w] = e1corr[i,w]/(1. - corr[i,w])
  e2corr[i,w] = e2corr[i,w]/(1. - corr[i,w])
  e1e1err[i,w] = e1e1err[i,w]/(1. - corr[i,w])
  e1e2err[i,w] = e1e2err[i,w]/(1. - corr[i,w])
  e2e2err[i,w] = e2e2err[i,w]/(1. - corr[i,w])


  ;; look for good g-band measurements

  w=where(struct.m_r[1] GT 0.0 AND struct.m_r[1] LT rcut, niband)

  w2=where(struct[w].m_e1e1err[1] LT errcut AND $
           struct[w].m_e2e2err[1] LT errcut, niband)
  w=w[w2]

  ;; correct g-band stuff
  e1corr[g,w]  = struct[w].m_e1_corr[1]
  e2corr[g,w]  = struct[w].m_e2_corr[1]
  e1e1err[g,w] = struct[w].m_e1e1err[1]
  e1e2err[g,w] = struct[w].m_e1e2err[1]
  e2e2err[g,w] = struct[w].m_e2e2err[1]

  e1corr[g,w] = e1corr[g,w]/(1. - corr[g,w])
  e2corr[g,w] = e2corr[g,w]/(1. - corr[g,w])
  e1e1err[g,w] = e1e1err[g,w]/(1. - corr[g,w])
  e1e2err[g,w] = e1e2err[g,w]/(1. - corr[g,w])
  e2e2err[g,w] = e2e2err[g,w]/(1. - corr[g,w])

END 

PRO test_cov3band,struct,w,$
                  cove,new_e1,new_e2,new_e1e1err,new_e1e2err,new_e2e2err,$
                  combine_flag,good

  IF n_elements(struct) EQ 0 THEN BEGIN 


      
      tags = ['m_e1','m_e2','m_e1e1err', 'm_e1e2err','m_e2e2err',$
              'm_e1_corr', 'm_e2_corr',$
              'm_rr_cc', 'm_rr_ccerr', $
              'm_cr4','m_r','petrocounts','reddening']

      run = 756
      rerun = 20
      FOR camcol=1,6 DO BEGIN 
          read_tsobj, [run,rerun,camcol], tstruct, /all, /corr, $
                      taglist=tags

          IF n_elements(struct) EQ 0 THEN BEGIN
              struct=temporary(tstruct) 
          ENDIF ELSE BEGIN
              concat_structs, temporary(struct),temporary(tstruct),tmp
              struct=temporary(tmp)
          ENDELSE 

      ENDFOR 
      struct.m_e1e1err = struct.m_e1e1err/2.
      struct.m_e1e2err = struct.m_e1e2err/2.
      struct.m_e2e2err = struct.m_e2e2err/2.

      return
  ENDIF 

  defval = -9999.
  rcut = 0.6
  errcut = 0.4
  minerr = 0.01

  cminrmag = 16.0
  cmaxrmag = 17.0
  print,'Magnitude range for corr matrix: ',cminrmag,cmaxrmag

  ;; pick objects in all bandpasses to measure
  ;; covariance matrix

  ;;
  rmag = struct.petrocounts[2]-struct.reddening[2]
  w=where(rmag LT cmaxrmag AND rmag GT cminrmag,nobj)

  ;; corrections applied, large objects
  w2=where(struct[w].m_r[1] GT 0.0 AND struct[w].m_r[1] LT rcut AND $
           struct[w].m_r[2] GT 0.0 AND struct[w].m_r[2] LT rcut AND $
           struct[w].m_r[3] GT 0.0 AND struct[w].m_r[3] LT rcut AND $
           struct[w].m_e2e2err[1] EQ struct[w].m_e2e2err[1] AND $
           struct[w].m_e2e2err[2] EQ struct[w].m_e2e2err[2] AND $
           struct[w].m_e2e2err[3] EQ struct[w].m_e2e2err[3], nobj)
  w=w[w2]

  w2=where(struct[w].m_e1e1err[1] LT errcut AND $
           struct[w].m_e2e2err[1] LT errcut AND $
           struct[w].m_e1e1err[2] LT errcut AND $
           struct[w].m_e2e2err[2] LT errcut AND $
           struct[w].m_e1e1err[3] LT errcut AND $
           struct[w].m_e2e2err[3] LT errcut, nobj)
  w=w[w2]

  clrs = [1,2,3]
  nclrs = n_elements(clrs)

  corr    = struct[w].m_r[clrs]
  e1corr  = struct[w].m_e1_corr[clrs]
  e2corr  = struct[w].m_e2_corr[clrs]
  e1e1err = struct[w].m_e1e1err[clrs]
  e1e2err = struct[w].m_e1e2err[clrs]
  e2e2err = struct[w].m_e2e2err[clrs]

  e1corr = e1corr/(1. - corr)
  e2corr = e2corr/(1. - corr)
  e1e1err = e1e1err/(1. - corr)
  e1e2err = e1e2err/(1. - corr)
  e2e2err = e2e2err/(1. - corr)


  ;; now just pick out good r-band galaxies: these will be our
  ;; base sample
  minrmag = 19.0
  maxrmag = 21.0
  print
  print,'Magnitude range for measurements: ',minrmag,maxrmag
  w=where(struct.m_r[2] GT 0.0 AND struct.m_r[2] LT rcut AND $
          struct.m_e1e1err[2] LT errcut AND $
          rmag LT maxrmag AND rmag GT minrmag, nrband)

  ;; now select particular objects (see program)
  test_cov3band_gete, struct[w], rcut, errcut, $
             e1corr, e2corr, e1e1err, e1e2err, e2e2err

  rsmear = struct[w].m_r[clrs]

  combine_ellip_cove1e2,e1corr,e2corr,e1e1err,e1e2err,e2e2err,rsmear,$
                        tnew_e1, tnew_e2, $
                        tnew_e1e1err, tnew_e1e2err, tnew_e2e2err,$
                        tcombine_flag, tgood,verbose=3,/besterr

  forprint,tnew_e1,new_e1,tnew_e2,new_e2,$
           tnew_e1e1err,new_e1e1err,tnew_e1e2err,new_e1e2err,$
           tnew_e2e2err,new_e2e2err


  ngood = n_elements(good)
  print,'Good measurements: '+ntostr(ngood)+'/'+ntostr(nrband)+' = '+$
        ntostr(float(ngood)/nrband)
  print,'Nbad = '+ntostr(nrband-ngood)+'/'+ntostr(nrband)+' = '+$
        ntostr(float(nrband-ngood)/nrband)

  g=0
  r=1
  i=2

  xrange=[-1,1]
  yrange=[-1,1]
  !p.multi=[0,0,2]
  plot,e1corr[r,good],new_e1[good]-e1corr[r,good],$
       psym=3,xrange=xrange,yrange=yrange
  oplot,[-1,1],[0,0],color=!blue

  plot,e2corr[r,good],new_e2[good]-e2corr[r,good],$
       psym=3,xrange=xrange,yrange=yrange
  oplot,[-1,1],[0,0],color=!blue

  print,mean(new_e1[good]-e1corr[r,good])
  print,mean(new_e2[good]-e2corr[r,good])

  !p.multi=0

END 
