;; http://www.ctio.noao.edu/~chrism/C4dr3/

FUNCTION c4::init, sample


  IF n_elements(sample) EQ 0 THEN BEGIN 
      message,'You must send the sample number',/inf
      message,"c4 = obj_new('c4', sample)"
  ENDIF 

  ;; The other defaults are mostly set in objshear::init.  See call below
  rmin = 20.0
  logbin = 1

  sample = fix(sample)

  type = 'c4'
  par = $
    {maskfile: '',      $
     sample: sample,    $
     catalog: '',       $
     source_sample:  1, $
     random_sample: -1, $
     max_allowed_angle: 6.0, $
     edgecuts: 1}


  CASE sample OF 

      1: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr3'
          par.random_sample = 9 ; The 10 Mpc sample
      END 
      2: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          par.catalog = 'dr3'
          par.random_sample = 9 ; The 10 Mpc sample
          par.max_allowed_angle = 20.0
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
  print,'max_allowed_angle =      ',par.max_allowed_angle
  print,'edgecuts =               ',par.edgecuts

  objshear_retval = $
    self->objshear::init(type, rmin, rmax, nbin, $
                         sigmacrit_style, shape_correction_style, $
                         logbin=logbin, par_struct=par)

  return, objshear_retval




  return,1
END 









FUNCTION c4::c4dir
  dir = expand_tilde(esheldon_config('lensinput_dir'))
  dir = concat_dir(dir,'c4')
  dir = concat_dir(dir,'catalog')
  return,dir
END 
FUNCTION c4::c4file
  dir=self->c4dir()
  file = 'sdss_c4_dr3_unpublished.fit'
  file = concat_dir(dir, file)
  return,file
END 
FUNCTION c4::get, rows=rows, columns=columns
  file=self->c4file()
  struct = mrdfits(file, 1, rows=rows, columns=columns)

  struct = rename_tags(struct, ['ra_mean', 'dec_mean'],['ra','dec'])
  add_tags, $
    struct, $
    ['clambda','ceta','loglum_r'],['0d','0d','-9999d'], $
    newstruct

  eq2csurvey, newstruct.ra, newstruct.dec, clam, ceta
  newstruct.clambda = clam
  newstruct.ceta = ceta

  w = where(newstruct.lum_r GT 0.0, nw)
  IF nw NE 0 THEN newstruct[w].loglum_r = alog10( newstruct[w].lum_r )
  return,newstruct
END 



PRO c4::loglum_bins, nbin, lowlim, highlim

  CASE nbin OF 
      2: BEGIN 
          ;; The 2/3 point in sorted lum_r
          lowlim  = [ 10.5,      11.494579 ]
          highlim = [ 11.494579, 12.27 ]
      END 
      ELSE: message,'Unsupported nbin: '+ntostr(nbin)
  ENDCASE 

END 
FUNCTION c4::loglum_where_string, nbin, labels=labels

  self->loglum_bins, nbin, lowlim, highlim

  delvarx, labels
  FOR i=0L, nbin-1 DO BEGIN 

      tstring = $
        '(struct.loglum_r GE '+ntostr(lowlim[i])+') AND '+$
        '(struct.loglum_r LT '+ntostr(highlim[i])+')'
      add_arrval,  tstring, where_string

      lowstr = ntostr(lowlim[i],5,/round)
      highstr = ntostr(highlim[i],5,/round)

      tlabel = lowstr+' < log(L!Dr!N) < '+highstr
      add_arrval, tlabel, labels

  END 
  return, where_string

END 


FUNCTION c4::subtype_nbin, subtype
  IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: nbin=c4->nbin(subtype)'
  ENDIF 
  CASE subtype OF
      'loglum2': return,2
      ELSE: message,'Unknown subtyep: '+strn(subtype)
  ENDCASE 
END 


FUNCTION c4::where_string, subtype, nbin=nbin, labels=labels
  IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: ws=c4->where_string(subtype)'
  ENDIF 

  nbin = self->subtype_nbin(subtype)

  CASE subtype OF 
      'loglum2': BEGIN 
          return, self->loglum_where_string(nbin, labels=labels)
      END 
      ELSE: message,'Unknown type: '+ntostr(subtype)
  ENDCASE 

END 

FUNCTION c4::labels, subtype, bin=bin
  t=self->where_string(subtype, labels=labels)
  IF n_elements(labels) EQ 0 THEN message,'labels undefined for subtype: '+$
    strlowcase(subtype)

  IF n_elements(bin) NE 0 THEN BEGIN 
      minb=min(bin, max=maxb)
      nl = n_elements(labels)
      IF minb LT 0 OR minb GE nl THEN message,'bin keyword out of range'
      labels = labels[bin]
  ENDIF 
  return,labels
END 










PRO c4::sub, subtype, randnum=randnum, zbinsize=zbinsize


  IF n_elements(subtype) EQ 0 THEN BEGIN 
      print,'-Syntax: c4->sub, subtype, randnum=, zbinsize='
      print
      message,'Halting'
  ENDIF 

  average_tags = ['z', $
                  'lum_r', 'loglum_r', $
                  'biwt1000', 'vmass1000']
  
  ;; just to make sure we know about this subtype
  nbin = self->subtype_nbin(subtype)

  self->objshear::sub_sample, subtype, randnum=randnum, $
    average_tags=average_tags, zhist=zhist, zbinsize=zbinsize

END 

PRO c4::plot_profile, type, subtype=subtype, bin=bin, $
                xlog=xlog, xmin=xmin, xmax=xmax, $
                ylog=ylog, ymin=ymin, ymax=ymax, $
                incolor=incolor, dops=dops
  
  IF n_elements(subtype) NE 0 THEN BEGIN 
      CASE subtype OF
          'loglum2': BEGIN 
              xlog=1
              ylog=1
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
              nsplit=2
              ymin = 0.1
              ymax = 900
          END 
          ELSE: message,'subtype: '+subtype+' not yet supported'
      ENDCASE  
  ENDIF ELSE BEGIN 
      ;; no sub sample
      CASE self->sample() OF 
          1: BEGIN 
              ymin = 0.1
              ymax = 100
          END 
          3: BEGIN 
              ymin = 1
              ymax = 100
          END 
          4: BEGIN 
              ymin = 1
              ymax = 100
          END 
          ELSE: BEGIN 
;              ymin = 1
;              ymax = 100
          END 
      ENDCASE 
      IF n_elements(xlog) EQ 0 THEN xlog=1
      IF n_elements(ylog) EQ 0 THEN ylog=1
  ENDELSE 


  self->objshear::plot_profile, type, subtype=subtype, bin=bin, $
    xlog=xlog, xmin=xmin, xmax=xmax, $
    ylog=ylog, ymin=ymin, ymax=ymax, $
    dops=dops, $
    charsize=charsize, label_charsize=label_charsize, $
    landscape=landscape, encapsulated=encapsulated, $
    incolor=incolor, $
    $
    nsplit=nsplit

END 




FUNCTION c4::cleanup
  return,1
END 
PRO c4__define

  struct = {$
             c4, $
             sample: 0, $
             INHERITS objshear $
           }

END 
