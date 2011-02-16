

PRO test_petrorad, str

  IF n_elements(str) EQ 0 THEN BEGIN 

      make_default_tags, tl

      taglist = ['ID',              $
                 'TYPE',            $
                 'PETROCOUNTS',     $
                 'PETROCOUNTSERR',  $
                 'STAR_L',          $
                 'EXP_L',           $
                 'DEV_L',           $
                 'PETRORAD',        $
                 'PETRORADERR',     $
                 'R_EXP',           $
                 'R_DEV']

      run=756
      camcol=3
      rerun=1
      start=100
      nf=100
      fetch_dir, run, camcol, rerun, dir

      read_tsobj, dir, str, start=start, nf=nf, taglist=taglist

  ENDIF 

  type=3
  maxpet=19.
  clr=2
  wexp = where(str.type[1] EQ type AND str.type[2] EQ type AND $
               str.exp_l[clr] GT str.dev_l[clr] AND $
               str.petrocounts[2] LT maxpet AND str.petrocounts[2] GT 0)
;  wexp = where(str.petrocounts[2] LT maxpet AND str.petrocounts[2] GT 0)


  wdev = where(str.type[1] EQ type AND str.type[2] EQ type AND $
               str.dev_l[clr] GT str.exp_l[clr] AND $
               str.petrocounts[2] LT maxpet AND str.petrocounts[2] GT 0)

      
  simpctable
  !p.background=!white
  !p.color=!black
  !p.multi=[0,0,2]

  sigma_clip,str[wexp].petrorad[clr]/str[wexp].r_exp[clr],mexp,sigexp,niter=3,nsig=3.5,/silent
  sigma_clip,str[wdev].petrorad[clr]/str[wdev].r_dev[clr],mdev,sigdev,niter=3,nsig=3.5,/silent

  print,mexp,mdev

  plot,str[wexp].petrorad[clr]/str[wexp].r_exp[clr],yrange=[0,10],psym=3
  oplot,[0,10000],[mexp,mexp]
  oplot,[0,10000],[4.43,4.43],line=2,color=!red
  plot,str[wdev].petrorad[clr]/str[wdev].r_dev[clr],yrange=[0,10],psym=3
  oplot,[0,10000],[mdev,mdev]
  oplot,[0,10000],[3.44,3.44],line=2,color=!red

  return
END 
