PRO test_myadmom, wh

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: test_myadmom, wh'
      return
  ENDIF 

  ;; ad_momi, im, incenx, inceny, shiftmax, sky, skysig, 
  ;;           ixx, iyy, ixy, rho4, uncert, wcenx, wceny, numiter, whyflag

  color_index=2
  
  setup_mystuff
;  bindir='/sdss3/usrdevel/philf/idl.lib/'
  bindir = !mybindir
  run=752
  rerun=1
  camcol=1
  field = 20
  nf=1

  maxmag = 22.
  ixx = fltarr(1)
  iyy = ixx
  ixy = ixx
  sky = ixx
  incenx = ixx
  inceny = ixx
  mag=ixx
  sky=ixx
  skysig=ixx
  numiter=lonarr(1)
  wcenx=ixx
  wceny=ixx
  whyflag=lonarr(1)
  rho4=ixx
  uncert=ixx
  n=long(1)


  fetch_dir, run, camcol, rerun, dir, atldir
  read_tsobj, dir, pstruct, start=field, nf=nf
  
  make_flag_struct,fs
  fs.satur='N'
  fs.satur_center='N'
  fs.bright = 'N'
  flag_select,pstruct,fs,color_index,_ss
  if(_ss[0] eq -1)then BEGIN
      print,'nogood ones'
      return
  ENDIF 
  pstruct=temporary( pstruct(_ss) )
  
  make_flag_struct, fs
  fs.blended='N'
  fs.moved='N'
  flag_select,pstruct,fs,color_index,_ss,/objc
  if(_ss[0] eq -1)then BEGIN
      print,'nogood ones'
      return
  ENDIF 
  pstruct=temporary( pstruct(_ss) )
  wt=where( (pstruct.petrocounts[2] -pstruct.reddening[2])LT maxmag,ntot)
;  ntot = n_elements(pstruct)
  pstruct=temporary(pstruct[wt])
  nobj = 300 < ntot

  print,'Hey ',ntot

  atlas_name, run, camcol, field, atlname
  atlname = atldir+atlname
  time = systime(1)

  IF wh EQ 1 THEN print, 'Measuring with phils!' ELSE print,'Measuring with mine!'
  print,'     iyy(fortran)      whyflag   iyy(IDL)      whyflag'
  philtime=fltarr(nobj)
  mytime = fltarr(nobj)
  ti=0L
  FOR kkk=4, nobj-1 DO BEGIN 

      im = 0
;      get_atlas_min, pstruct[kkk].id, color_index, atlname, $
;        im, s, row0, col0

      setzero,im
      get_atlas,pstruct,kkk,dir=atldir,clr=2,imr=im,/noprompt,/nodisplay,$
        col0=col0,row0=row0,/silent

      IF 1 THEN BEGIN 

          im = float(im)
          sz = size(im)
          nx = sz[1]
          ny = sz[2]

          ixx[0] = 4.5 & tixx=4.5
          iyy[0] = 4.5 & tiyy=4.5
          ixy[0] = 0.0 & tixy=0.0
          sky[0] = 1000.
          skysig[0]=5.

          shiftmax = 2.*pstruct[kkk].petrorad(color_index) > 2.
          shiftmax = shiftmax < 10.

          incenx[0] = pstruct[kkk].colc(color_index)-col0+0.5
          inceny[0] = pstruct[kkk].rowc(color_index)-row0+0.5
          twhyflag=0

;          IF wh EQ 1 THEN BEGIN 
              time=systime(1)
              ff=call_external(bindir+'ad_momi.so','ad_mom_', $
                               incenx, inceny, $
                               ixx, iyy, ixy, $
                               n, shiftmax, $
                               im, sky, nx, ny, uncert, skysig, mag,$
                               numiter, wcenx, wceny, whyflag, rho4)
              philtime[ti] = systime(1)-time
;          ENDIF ELSE BEGIN
              time=systime(1)
              ad_momi2, im, incenx, inceny, shiftmax, $
                       sky, skysig, $
                       tixx, tiyy, tixy, rho4, tuncert, $
                       wcenx, wceny, numiter, twhyflag
              mytime[ti] = systime(1)-time
;          ENDELSE 
          ti=ti+1L
      ENDIF ELSE print,'nogood'

;      IF ixx[0] NE 0. THEN print,(ixx[0]-tixx)/ixx[0] ELSE print,'ixx=0'
      print,iyy[0],whyflag[0],tiyy,twhyflag
;      IF (whyflag[0] NE 0) OR (twhyflag[0] NE 0) THEN print,iyy[0],whyflag[0],tiyy,twhyflag
;      print,iyy[0],whyflag[0],tiyy,twhyflag

  ENDFOR 
  ptime,total(philtime)
  ptime,total(mytime)


  return
END 
