FUNCTION bax_lensing::init, sample

  IF n_elements(sample) EQ 0 THEN BEGIN 
      message,'You must send the sample number',/inf
      message,"ll = obj_new('lrg_lensing', sample)",/inf
      message,''
  ENDIF 

  rmin = 20.0
  logbin = 1

  sample = fix(sample)

  type = 'bax'
  par = $
    { $
      maskfile: '', $
      sample: sample, $
      source_sample: 1, $
      random_sample: -1, $
      catalog: '',$
      zbuffer: 0.0 $
    }

  CASE sample OF 
      1: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'bax'
          par.random_sample = 9 ; The 10 Mpc sample
          par.source_sample = 1 ; old cuts just for a check
      END 
      2: BEGIN 
          rmax = 36567.0
          nbin = 21
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'bax'
          par.random_sample = 10 ; The 30 Mpc sample
          par.source_sample = 1 ; old cuts just for a check
      END 
      ELSE: message,'Unknown sample: '+ntostr(sample)
  ENDCASE 

  self.sample = sample

  print,'type   =                 ',type
  print,'sample =                 ',sample,format='(A,I7)'
  print,'source_sample =          ',par.source_sample,format='(A,I7)'
  print,'random_sample =          ',par.random_sample,format='(A,I7)'
  print,'rmin =                   ',rmin
  print,'rmax =                   ',rmax
  print,'nbin =                   ',nbin
  print,'sigmacrit_style =        ',sigmacrit_style,format='(A,I7)'
  print,'shape_correction_style = ',shape_correction_style,format='(A,I7)'
  objshear_retval = $
    self->objshear::init(type, rmin, rmax, nbin, $
                         sigmacrit_style, shape_correction_style, $
                         logbin=logbin, par_struct=par)

  return, objshear_retval



END 

FUNCTION bax_lensing::catdir
  dir = esheldon_config('lensinput_dir')
  dir = concat_dir(dir, 'bax/catalog')
  return,dir
END 

FUNCTION bax_lensing::catname, original=original
  dir = self->catdir()
  file = concat_dir(dir, 'bax.fit')

  IF NOT keyword_set(original) THEN BEGIN  
      file = repstr(file, '.fit', '_fix.fit')
  ENDIF 
  return,file
END 

FUNCTION bax_lensing::get, original=original
  file = self->catname(original=original)
  print
  print,'Reading file: ',file
  struct = mrdfits(file, 1)
  return,struct
END 











FUNCTION bax_lensing::zcuts

  ;; use with strict gt and lt
  return,[0.05156, 0.3]

END 

FUNCTION bax_lensing::zselect, z, nkeep

  zcuts = self->zcuts()
  keep  = where(z GT zcuts[0] AND z LT zcuts[1], nkeep)
  return, keep

END 



PRO bax_lensing::fix

  outfile = self->catname()
  print
  print,'Will write to new file: ',outfile
  
  print

  orig = self->get(/original)
  norig = n_elements(orig)

  ;; Old catalogs didn't have clambda,ceta
  IF NOT tag_exist(orig[0], 'clambda') THEN BEGIN 
      print
      print,'Adding clambda,ceta'
      newstruct = create_struct(orig[0], 'clambda', 0d, 'ceta', 0d)
      newstruct = replicate(newstruct, norig)

      copy_struct, orig, newstruct

      eq2csurvey, newstruct.ra, newstruct.dec, clambda, ceta
      newstruct.clambda = clambda
      newstruct.ceta = ceta
  ENDIF ELSE BEGIN 
      newstruct = temporary(orig)
  ENDELSE 

  print
  print,'Applying redshift cuts'
  keep2 = self->zselect(newstruct.z, nkeep2)
  newstruct = newstruct[keep2]
  print,'Kept '+ntostr(nkeep2)+'/'+ntostr(norig)



  print
  print,'Writing new file: '+outfile
  mwrfits, newstruct, outfile, /create


END 

FUNCTION bax_lensing::subtype_nbin, subtype

  CASE subtype OF 
      'tx2': return,2
      ELSE: message,'Unknown subtype: '+strn(subtype)
  ENDCASE 

END 

PRO bax_lensing::tx_bins, nbin, lowlim, highlim

  CASE nbin OF 
      2: BEGIN 
          lowlim = [0.9, 5.37]
          highlim = [5.37, 9.6]
      END 
      ELSE: message,'Unknown nbin: '+strn(nbin)
  ENDCASE 

END 

FUNCTION bax_lensing::tx_where_string, nbin

  self->tx_bins, nbin, lowlim, highlim

  IF nbin EQ 1 THEN where_string = '' $
  ELSE where_string = strarr(nbin)
  
  FOR i=0L, nbin-1 DO BEGIN 

      where_string[i] = $
        'struct.tx GE '+ntostr(lowlim[i])+' AND '+$
        'struct.tx LT '+ntostr(highlim[i])+' AND '+$
        'struct.lx gt 0'

  ENDFOR 
  return,where_string

END 


FUNCTION bax_lensing::where_string, subtype, nbin=nbin
  IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: ws=bax->where_string(subtype)'
  ENDIF 

  nbin = self->subtype_nbin(subtype)

  CASE subtype OF 
      'tx2': BEGIN 
          return, self->tx_where_string(nbin)
      END 
      ELSE: message,'Unknown type: '+ntostr(subtype)
  ENDCASE 

END 


PRO bax_lensing::sub, subtype, randnum=randnum, zbinsize=zbinsize

  average_tags = ['z', 'lx', 'tx']
  
  zhist = 1
  IF n_elements(zbinsize) EQ 0 THEN zbinsize = 0.03
  self->objshear::sub_sample, subtype, randnum=randnum, $
    average_tags=average_tags, zhist=zhist, zbinsize=zbinsize

END 









PRO bax_lensing::plot, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      plotdir = self->plot_dir(/createdir)

      file = self->type()+'_deltasig_'+self->sample_string()+'.eps'
      file = concat_dir(plotdir, file)

      begplot, file, /color, /encapsulated
  ENDIF 

  t = self->lensread('corr')

  plot_density_contrast, t, /log, /mpc, aspect=!gratio

  IF keyword_set(dops) THEN endplot, /trim_bbox

END 


PRO bax_lensing::tx_plot, subtype, cumulative=cumulative, dops=dops

  IF n_params() LT 1 THEN BEGIN 
      message,'-syntax: bax->tx_plot, subtype, /cumulative, /dops'
  ENDIF 

  IF keyword_set(dops) THEN BEGIN 
      plotdir = self->plot_dir(subtype=subtype, /createdir)

      file = 'bax_'+self->sample_string()+'_'+subtype

      IF keyword_set(cumulative) THEN BEGIN 
          file = file + '_cumulative'
      ENDIF 
      file = file + '.eps'

      file = concat_dir(plotdir, file)

      begplot, file, /color, /encapsulated
  ENDIF 

  esheldon_setup

  sh =self->lensread('corr')
  tx = self->lensread('corr', subtype=subtype)

  clr0 = !darkGreen
  clr1 = !red

 
  ytitle = !deltaytitle
  xtitle = !mpcxtitle2

  IF keyword_set(cumulative) THEN BEGIN 
      shtag = where(tag_names(sh) EQ 'TSIGMA')
      shetag = where(tag_names(sh) EQ 'TSIGMAERR')
      txtag = where(tag_names(tx) EQ 'TSIGMA')
      txetag = where(tag_names(tx) EQ 'TSIGMAERR')
      
      addstr = 'Cumulative '
  ENDIF ELSE BEGIN 
      shtag = where(tag_names(sh) EQ 'SIGMA')
      shetag = where(tag_names(sh) EQ 'SIGMAERR')
      txtag = where(tag_names(tx) EQ 'SIGMA')
      txetag = where(tag_names(tx) EQ 'SIGMAERR')
      addstr = ''
  ENDELSE 
  ytitle = addstr + ytitle
  erase & multiplot, [1,2]

  ploterror, sh.meanr/1000, sh.(shtag), sh.(shetag), $
    /xlog, /ylog,yrange=[1,1.e3],ystyle=3,xrange=[0.02,20],xstyle=3,$
    ytitle=ytitle, ytickf='loglabels'
  oploterror, tx[1].meanr/1000, tx[1].(txtag), tx[1].(txetag), $
    color=clr1
  oploterror, tx[0].meanr/1000, tx[0].(txtag), tx[0].(txetag), $
    color=clr0

  units = '5.37 keV'
  legend, $
    ['all', 'L!DX!N < '+units, 'L!DX!N > '+units], $
    line=0, color=[!p.color, clr0, clr1], /right, box=0, charsize=1

  multiplot

  ratio0 = tx[0].tsigma/sh.tsigma
  ratio0err = ratio0*sqrt( (sh.(shetag)/sh.(shtag))^2 + $
                           (tx[0].(txetag)/tx[0].(txtag))^2 )

  ratio1 = tx[1].tsigma/sh.tsigma
  ratio1err = ratio1*sqrt( (sh.(shetag)/sh.(shtag))^2 + $
                           (tx[1].(txetag)/tx[1].(txtag))^2 )


  plot, [0], /nodata, $
    /xlog, yrange=[0,2],ystyle=3,xrange=[0.02,20],xstyle=3, $
    xtitle=xtitle, ytitle='L!DX!N!Ubin!N / all', $
    xtickf='loglabels'

  oplot, [0.001, 10000], [1,1]
  oploterror, sh.meanr/1000, ratio0, ratio0err, $
    color=clr0
  oploterror, sh.meanr/1000, ratio1, ratio1err, $
    color=clr1
    


  multiplot, /reset

  IF keyword_set(dops) THEN endplot, /trim_bbox


END 







PRO bax_lensing__define

  struct = { $
             bax_lensing, $
             sample: 0, $
             INHERITS objshear $
           }

END 
