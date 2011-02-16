PRO map_ellip, e1, e2, uncert, mag, klist=klist, elist=elist, nosex=nosex

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: map_ellip, e1, e2, uncert, mag, klist=klist, elist=elist, nosex=nosex'
      return
  ENDIF 
  
  IF n_elements(klist) EQ 0 THEN BEGIN 

      start = 5
      max   = 23
      nlist = max - start + 1

      dir = '/sdss4/data1/esheldon/GAL_GAL/MASSMAPS/'
      base = dir+'galaxies_S10_R100_'
      
      klist = strarr(nlist)
      elist = klist

      ii = ntostr( start+indgen(max-start+1) )
      klist = base + 'kappa_N'+ii+'.fit'
      elist = base + 'kerr_N'+ii+'.fit'

  ENDIF 
  
  help,klist
  nlist =  n_elements(klist)
  IF NOT keyword_set(nosex) THEN BEGIN 
      
      e1 = fltarr(nlist)
      e2 = e1
      e  = e1
      uncert = e1
      mag = e1
      FOR i=0, nlist-1 DO BEGIN 

          im = mrdfits(klist[i], /silent)
          err = mrdfits(elist[i], /silent)

          imsig = im/err

          sdss_extract, imsig, tmp

          mag[i] = max(tmp.mag_best)

          w=where(tmp.mag_best EQ mag[i])

          e1[i] = tmp[w].e1_ad
          e2[i] = tmp[w].e2_ad
          uncert[i] = tmp[w].uncert_ad

      ENDFOR 

      mag = 30. - mag
  ENDIF 
  
  e = sqrt( e1^2 + e2^2 )
  w = where( e1 NE -10.,  nw)
  nbad = nlist - nw
  IF nbad NE 0 THEN BEGIN 
      print,'Found ',ntostr(nbad),' not converged'
  ENDIF 

  IF nw NE 0 THEN BEGIN 
      
      wmom, e[w],  sqrt(2.)*uncert[w], wme, wsige, werre
      wmom, e1[w], uncert[w], wme1, wsige1, werre1
      wmom, e2[w], uncert[w], wme2, wsige2, werre2

      mess_e  = ['Mean e:  ' + ntostr(wme)  + ' +/- ' + ntostr(werre), $
                 'Stdv e:  ' + ntostr(wsige) ]
      mess_e1 = ['Mean e1: ' + ntostr(wme1) + ' +/- ' + ntostr(werre1), $
                 'Stdv e1:  ' + ntostr(wsige1) ]
      mess_e2 = ['Mean e2: ' + ntostr(wme2) + ' +/- ' + ntostr(werre2), $
                 'Stdv e2:  ' + ntostr(wsige2) ]

      ytit = 'Number'
      xtit = 'e'
      plothist, e[w], bin=.05, ytit=ytit, xtit=xtit
      legend, mess_e
      key=get_kbrd(1)

      xtit = 'e1'
      plothist, e1[w], bin=.05, ytit=ytit, xtit=xtit
      legend, mess_e1
      key=get_kbrd(1)
  
      xtit = 'e2'
      plothist, e2[w], bin=.05, ytit=ytit, xtit=xtit
      legend, mess_e2
      key=get_kbrd(1)

      xtit = 'uncert'
      plothist, uncert[w], bin=.02, ytit=ytit, xtit=xtit
      key=get_kbrd(1)

      xtit = 'mag_best'
      plothist, mag[w], bin=.1, ytit=ytit, xtit=xtit


  ENDIF 

return
END 
