PRO get_maps, extension, clr, k, kerr, ksig, n, nerr, nsig, kcat=kcat, ncat=ncat, step=step, bin=bin, slength=slength, med=med, maxx=maxx, minx=minx, maxy=maxy, miny=miny, name=name

  IF N_params() LT 2 THEN BEGIN 
     print,'-Syntax: get_maps, extension, clr, k, kerr, ksig, n, nerr, nsig, step=step, bin=bin, S=S, kcat=kcat, ncat=ncat,maxx=maxx, minx=minx, maxy=maxy, miny=miny, name=name'
     print,''
     print,'Use doc_library,"get_maps"  for more help.'  
     return
  ENDIF 

  dir='/sdss4/data1/esheldon/LARGEMAPS/'
  ext = ntostr(long(extension))
  colors = ['u','g','r','i','z']

  k=mrdfits(dir+'big_map_kappa_'+colors[clr]+'_N'+ext+'.fit')
  kerr=mrdfits(dir+'big_map_kerr_'+colors[clr]+'_N'+ext+'.fit')
    
  n=mrdfits(dir+'big_map_noise_'+colors[clr]+'_N'+ext+'.fit')
  nerr=mrdfits(dir+'big_map_nerr_'+colors[clr]+'_N'+ext+'.fit')
    
  ksig = k/kerr
  nsig = n/nerr

  nstep = n_elements(step)
  IF nstep NE 0 THEN BEGIN
      print,'Doing Analysis'

      rangefile = dir+'big_map_range_'+colors[clr]+'_N'+ext+'.txt'
      readcol, rangefile, char1, char2, v, format='A4, A1, F'
      minx = v[0]
      maxx = v[1]
      miny = v[2]
      maxy = v[3]

      sz=size(k)
      sx=sz[1]
      sy=sz[2]
      nstep = sy/sx/step

      ;; Only redo if cat's don't exist
      nkcat = n_elements(kcat) & nncat = n_elements(ncat)

      IF (nkcat EQ 0) OR (nncat EQ 0) THEN BEGIN 
          mextract, k, ksig, kcat, n, nsig, ncat, step=step
      ENDIF 
      
      IF n_elements(name) NE 0 THEN begplot,name=name+'.ps'

      masshist, kcat, ncat, bin, maxx=maxx, minx=minx, maxy=maxy, miny=miny, $
                step=step, nstep=nstep, slength=slength, name=name, med=med

      IF n_elements(name) NE 0 THEN ep
      
  ENDIF 

return
END 

