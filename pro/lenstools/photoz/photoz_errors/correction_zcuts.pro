PRO correction_zcuts

  dir = '/net/cheops1/data0/esheldon/test_sigmacriterr/inputz/'
  cd,dir

  xrange = [0.0, 0.3]
  ;;yrange = [8.0, 12.0]
  xtitle = 'z!DL!N'
  ytitle = 'Correction Factor - 1'

  tend = ['N1.fit', 'N2.fit']
  type = ['1 '+!csym.sigma, '2 '+!csym.sigma]
  psyms = [1, 4]
  colors = [!p.color, !red]
  n_end = n_elements(tend)

  FOR ie=0L, n_end-1 DO BEGIN 

      delvarx, zL, mean_mdenscont, sig_mdenscont, $
               err_mdenscont, mean_mdensconttrue

      files = findfile("fits_*"+tend[ie])
            
      nf = n_elements(files)
      
      FOR i=0L, nf-1 DO BEGIN 
          
          t=mrdfits(files[i],1)
          
          add_arrval, t.zL[0], zL
          add_arrval, t.mean_mdenscont, mean_mdenscont
          add_arrval, t.sig_mdenscont, sig_mdenscont
          add_arrval, t.err_mdenscont, err_mdenscont
          add_arrval, t.mean_mdensconttrue, mean_mdensconttrue
          
      ENDFOR 
      
      rat = mean_mdensconttrue/mean_mdenscont
      raterr = rat*err_mdenscont/mean_mdenscont
      rat = rat-1

      IF ie EQ 0 THEN BEGIN 
xrange=[0.04,0.3]
yrange=[1.005,1.2] - 1.
          aploterror, 1, zL, rat, raterr, psym=psyms[ie], $
                      yrange=yrange, xrange=xrange, xstyle=1,ystyle=1, $
                      xtitle=xtitle, ytitle=ytitle, $
                      /xlog,/ylog
          add_labels,xtickv=[0.05,0.2]

          ;;oplot, zL, rat, line=0
          ;;oplot, zl, mean_mdensconttrue

          fitpower, zL, rat, replicate(1.0,nf), [0.1,1.0], yfit, a
          oplot, zL, yfit, color=colors[ie]

      ENDIF ELSE BEGIN 
          oploterror, zL, rat, raterr, psym=psyms[ie], $
                      color=colors[ie],errcolor=colors[ie], line=0
          ;;oplot, zL, rat, line=0, color=colors[ie]

          fitpower, zL, rat, replicate(1.0,nf), [0.1,1.0], yfit, a
          oplot, zL, yfit, color=colors[ie]

;          fitlin, zL, rat, raterr, a, siga, b, sigb
;          a = linfit(zL, rat)
;          print,a
          ;;oplot, zL, a + b*zL
;          oplot, zL, a[0] + a[1]*zL

      ENDELSE 

  ENDFOR 

  legend, type, psym=psyms, color=colors, thick=replicate(!p.thick,n_end)

END 
