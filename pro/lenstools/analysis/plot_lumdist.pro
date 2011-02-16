PRO plot_lumdist, usemultiplot=usemultiplot, linefills=linefills

  massclr=2
  basedir =  '/sdss5/data0/lensout/'
  indir = basedir+'luminosities/'
  infile = indir + 'lumdist_data_'+!colors[massclr]+'lensing.fit'

  lumd=mrdfits(infile, 1)

  bin = 0.1
  ticklen=0.04

  xrange=fltarr(5,2)
  xrange[0,*] = [-22,-14]
  xrange[1,*] = [-24,-16]
  xrange[2,*] = [-24,-14]
  xrange[3,*] = [-25,-16]
  xrange[4,*] = [-25,-16]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Use multiplot
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(usemultiplot) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Use line filling
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF keyword_set(linefills) THEN BEGIN 

          outfile = indir + 'lumdist_multiplot_linefill_'$
            +!colors[massclr]+'lensing.ps'
          begplot,name=outfile

          erase & multiplot,[2,3]
          FOR clr=0L, 4 DO BEGIN 
              xtitle='' & ytitle=''
              IF clr EQ 4 THEN BEGIN
                  xtitle='M - 5 log!D10!N h'
                  ytitle = 'Number of Galaxies'
              ENDIF 
              ;;IF clr EQ 3 THEN BEGIN
              ;;    xtitle='M - 5 log!D10!N h'
              ;;ENDIF 
              xrange=[-25,-14]
              yrange=[0, 1700]
              cs=!colors[clr]
              command='ab1=lumd.'+cs+'absmag1 & ab2=lumd.'+cs+'absmag2 & ab3=lumd.'+cs+'absmag3 & '+$
                'ab4=lumd.'+cs+'absmag4'
              IF NOT execute(command) THEN message,'Error'
              
              thick=3
              plothist, ab1, bin=bin, xrange=xrange, yrange=yrange,$
                xtitle=xtitle,$
                ytitle=ytitle, thick=thick, $
                /fill,/fline,forientation=135,linestyle=1,$
                ystyle=1, noclip=0,xstyle=1, ticklen=ticklen
              plothist, ab2, bin=bin, /overplot,$
                /fill,/fline,forientation=45, thick=thick, noclip=0
              plothist, ab3, bin=bin, /overplot,$
                /fill,/fline,forientation=0, thick=thick, noclip=0
              plothist, ab4, bin=bin, /overplot,$
                /fill,/fline,forientation=135, thick=thick, noclip=0
              legend,!colorsp[clr],/right,box=0
              
              doxaxis=0
              IF clr EQ 2 OR clr EQ 3 THEN doxaxis=1
              IF clr NE 4 THEN multiplot,doxaxis=doxaxis
          ENDFOR 
          
          multiplot,/reset
      ENDIF ELSE BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; No line filling
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          outfile = indir + 'lumdist_multiplot_'$
            +!colors[massclr]+'lensing.ps'
          begplot,name=outfile

          erase & multiplot,[2,3]
          FOR clr=0L, 4 DO BEGIN 
              xtitle='' & ytitle=''
              IF clr EQ 4 THEN BEGIN
                  xtitle='M - 5 log!D10!N h'
                  ytitle = 'Number of Galaxies'
              ENDIF 
              ;;IF clr EQ 3 THEN BEGIN
              ;;    xtitle='M - 5 log!D10!N h'
              ;;ENDIF 
              xrange=[-25,-14]
              yrange=[0, 1700]
              cs=!colors[clr]
              command='ab1=lumd.'+cs+'absmag1 & ab2=lumd.'+cs+'absmag2 & ab3=lumd.'+cs+'absmag3 & '+$
                'ab4=lumd.'+cs+'absmag4'
              IF NOT execute(command) THEN message,'Error'
              
              thick=3
              plothist, ab1, bin=bin, xrange=xrange, yrange=yrange,$
                xtitle=xtitle,$
                ytitle=ytitle, thick=thick, $
                linestyle=0,$
                ystyle=1, xstyle=1, ticklen=ticklen
              plothist, ab2, bin=bin, /overplot,$
                linestyle=2, thick=thick
              plothist, ab3, bin=bin, /overplot,$
                linestyle=3, thick=thick
              plothist, ab4, bin=bin, /overplot,$
                linestyle=4, thick=thick
              legend,!colorsp[clr],/right,box=0
              
              doxaxis=0
              IF clr EQ 2 OR clr EQ 3 THEN doxaxis=1
              IF clr NE 4 THEN multiplot,doxaxis=doxaxis
          ENDFOR 
          
          multiplot,/reset

      ENDELSE 
  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; No multiplot
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Use line filling
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF keyword_set(linefills) THEN BEGIN 

          outfile = indir + 'lumdist_linefill_'$
            +!colors[massclr]+'lensing.ps'
          begplot,name=outfile

          !p.multi=[0,2,3]
          FOR clr=0L, 4 DO BEGIN 

              ;;xtitle='M('+!colorsp[clr]+') - 5 log!D10!N h'
              xtitle='M - 5 log!D10!N h'
              ytitle = 'Number of Galaxies'

              cs=!colors[clr]
              command='ab1=lumd.'+cs+'absmag1 & ab2=lumd.'+cs+'absmag2 & ab3=lumd.'+cs+'absmag3 & '+$
                'ab4=lumd.'+cs+'absmag4'
              IF NOT execute(command) THEN message,'Error'
              
              thick=3
;              plothist, ab1, bin=bin,$
;                xtitle=xtitle,$
;                ytitle=ytitle, thick=thick, $
;                /fill,/fline,forientation=135,linestyle=1,$
;                /ylog,yrange=[1,10000],$
;                charsize=2,xrange=xrange[clr,*], noclip=0, ticklen=ticklen
              plothist, ab1, bin=bin,$
                xtitle=xtitle,$
                ytitle=ytitle, thick=thick, $
                /fill,/fline,forientation=135,linestyle=1,$
                charsize=2,xrange=xrange[clr,*], noclip=0, ticklen=ticklen

              plothist, ab2, bin=bin, /overplot,$
                /fill,/fline,forientation=45, thick=thick, noclip=0
              plothist, ab3, bin=bin, /overplot,$
                /fill,/fline,forientation=0, thick=thick, noclip=0
              plothist, ab4, bin=bin, /overplot,$
                /fill,/fline,forientation=135, thick=thick, noclip=0
              legend,!colorsp[clr],/right,box=0
              
          ENDFOR 
          
          !p.multi=0
      ENDIF ELSE BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; No line filling
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          outfile = indir + 'lumdist_'$
            +!colors[massclr]+'lensing.ps'
          begplot,name=outfile

          !p.multi=[0,2,3]
          FOR clr=0L, 4 DO BEGIN 

              ;;xtitle='M('+!colorsp[clr]+') - 5 log!D10!N h'
              xtitle='M - 5 log!D10!N h'
              ytitle = 'Number of Galaxies'

              cs=!colors[clr]
              command='ab1=lumd.'+cs+'absmag1 & ab2=lumd.'+cs+'absmag2 & ab3=lumd.'+cs+'absmag3 & '+$
                'ab4=lumd.'+cs+'absmag4'
              IF NOT execute(command) THEN message,'Error'
              
              thick=3
              plothist, ab1, bin=bin,$
                xtitle=xtitle,$
                ytitle=ytitle, thick=thick, $
                linestyle=0,$
                charsize=2,xrange=xrange[clr,*],/ylog, ticklen=ticklen
              plothist, ab2, bin=bin, /overplot,$
                linestyle=2, thick=thick
              plothist, ab3, bin=bin, /overplot,$
                linestyle=3, thick=thick
              plothist, ab4, bin=bin, /overplot,$
                linestyle=4, thick=thick
              legend,!colorsp[clr],/right,box=0
              
          ENDFOR 
          
          !p.multi=0

      ENDELSE 
      

  ENDELSE 
  endplot
;  !p.multi=0

END 
