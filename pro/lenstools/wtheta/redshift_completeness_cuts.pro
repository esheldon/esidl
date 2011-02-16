PRO redshift_completeness_cuts,clr

  nn=10000L
  z=arrscl(findgen(nn), 0.01, .4)
  a=angdist_lambda(z,dlum=dlum) ;Mpc
  dlum = dlum*1.e6              ;pc

  fraclstar = 0.1
  max_mag = [-1000., 21.0, 21.0, 21.0, 19.8]

  Mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  Msun = [6.39,5.07,4.62,4.52,4.48]
  Lstar = 10.0^((Mstar-Msun)/(-2.5))

;  print,lstar

  karray=mrdfits('/sdss3/data5/mckay/fukugita_k_corrections_array.fit',0)
  tkarray = reform( karray[clr,*,1:6] )
  ;;This array has three indices:
  ;;	First: Color u,g,r,i,z
  ;;	Second: Redshift bin
  ;;	Third: z, then six types (E,SO,Sab,Sbc,Scd,Im)
  ;;
  ;;	For example karray(2,*,0) is the redshift values for the r data
  ;;		    karray(2,*,2) is the r' k corrections for SO galaxies
  
  ;;First set up the color information required to type the galaxies
  gmr_color=karray(1,*,1:6)-karray(2,*,1:6)
  zval=reform( karray(1,*,0) )
  gmr_color=reform(gmr_color)

  ;; Ellipticals have largest k-correction
  type = 0L

  kcorr_ellip = -1.0*interpolate(tkarray, z/0.05, replicate(type,nn))

  m_max = Mstar[clr] - 2.5*alog10(fraclstar) + $
    5.0*alog10(dlum/10.0) - kcorr_ellip

  maxz = max( z[ where(m_max LE max_mag[clr]) ] )
  print,'For max_mag['+!colors[clr]+'] = '+ntostr(max_mag[clr])+' max Z = '+ntostr(maxz)
  plot, z, m_max,ytitle=!colors[clr]+'-band mag',xtitle='Z'
  oplot,[maxz, maxz], [-10, 1000], line=2
  oplot,[0,10000],[max_mag[clr],max_mag[clr]],line=2
  legend,'max_mag = '+ntostr(max_mag[clr]),/right,box=0

END   
