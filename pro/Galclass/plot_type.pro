PRO plot_type, plot=plot

;  IF n_params() LT 1 THEN BEGIN 

;  ENDIF 
  
  IF NOT keyword_set(plot) THEN plot=0

  indir = '/sdss4/data1/esheldon/GAL_GAL/spectra/'
  
  fg = indir + [$
                 'spiral_752_756_g_sum_N1.fit' $
                 ,'ellip_752_756_g_sum_N1.fit'  $
                 ]
  fr = indir + [$
                 'spiral_752_756_r_sum_N1.fit' $
                 ,'ellip_752_756_r_sum_N1.fit'  $
                 ]
  fi = indir + [$
                 'spiral_752_756_i_sum_N1.fit' $
                 ,'ellip_752_756_i_sum_N1.fit'  $
                 ]

  grand = indir + [$
                    'spiral_rand_752_756_g_sum_N1.fit' $
                    ,'ellip_rand_752_756_g_sum_N1.fit'  $
                 ]
  rrand = indir + [$
                    'spiral_rand_752_756_r_sum_N1.fit' $
                    ,'ellip_rand_752_756_r_sum_N1.fit'  $
                 ]
  irand = indir + [$
                    'spiral_rand_752_756_i_sum_N1.fit' $
                    ,'ellip_rand_752_756_i_sum_N1.fit'  $
                 ]

  gzf = indir + [$
                  'spiral_752_756_g_z_N1.fit' $
                  ,'ellip_752_756_g_z_N1.fit' $
                ]
  rzf = indir + [$
                  'spiral_752_756_r_z_N1.fit' $
                  ,'ellip_752_756_r_z_N1.fit' $
                ]
  izf = indir + [$
                  'spiral_752_756_i_z_N1.fit' $
                  ,'ellip_752_756_i_z_N1.fit' $
                ]

  glensumf = indir + [$
                  'spiral_752_756_g_lensum_N1.fit' $
                  ,'ellip_752_756_g_lensum_N1.fit' $
                ]
  rlensumf = indir + [$
                  'spiral_752_756_r_lensum_N1.fit' $
                  ,'ellip_752_756_r_lensum_N1.fit' $
                ]
  ilensumf = indir + [$
                  'spiral_752_756_i_lensum_N1.fit' $
                  ,'ellip_752_756_i_lensum_N1.fit' $
                ]

  nf=n_elements(fg)

  FOR i=0, nf-1 DO BEGIN 
      
      tmpstr = str_sep(fg[i], 'g_sum_N')
      tmpstr2 = str_sep(tmpstr[1], '.fit')
      name=tmpstr[0]+'plots_N'+tmpstr2[0]+'.ps'
      
      IF plot THEN begplot, name=name,/color

      IF i EQ 1 THEN BEGIN 
          sigvrange=[50., 250.]
          cutrange=[100., 1300.]
          normrange=[0., 15.]
          powrange=[0.3, 1.4]
      ENDIF ELSE BEGIN
          sigvrange=[0., 200.]
          cutrange=[100., 1300.]
          normrange=[0., 6.]
          powrange=[0.3, 1.4]
      ENDELSE 
      plotclrzshear, fg[i], fr[i], fi[i],$ 
        grand=grand[i],rrand=rrand[i],irand=irand[i],$
        gzfiles=gzf[i], rzfiles=rzf[i], izfiles=izf[i], $
        munit=1.e12, $
        sigvrange=sigvrange, cutrange=cutrange, $
        normrange=normrange, powrange=powrange;, $
        ;glensum=glensumf, rlensum=rlensumf, ilensum=ilensumf

      IF plot THEN ep
      IF NOT plot THEN BEGIN
          key=get_kbrd(1)
          IF key EQ 'q' THEN return
      ENDIF 
      
  ENDFOR 

return
END 
