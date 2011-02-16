FUNCTION vagc_lensing::init, sample

  IF n_elements(sample) EQ 0 THEN BEGIN 
      message,'You must send the sample number',/inf
      message,"vl = obj_new('vagc_lensing', sample)"
  ENDIF 

  ;; The other defaults are mostly set in objshear::init.  See call below
  rmin = 20.0
  logbin = 1

  sample = fix(sample)

  type = 'vagc'
  par = $
    {maskfile: '',      $
     sample: sample,    $
     catalog: '',       $
     source_sample:  1, $
     random_sample: -1, $
     max_allowed_angle: 6.0, $
     edgecuts: 1}


  CASE sample OF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; lssfull0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      1: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'lssfull0'
          par.random_sample = 9 ; The 10 Mpc sample
      END 
      2: BEGIN 
          ;; NO EDGE CUTS!! 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'lssfull0'
          par.random_sample = 9 ; The 10 Mpc sample
          par.max_allowed_angle = 20.0
          par.edgecuts = 0

      END 

      3: BEGIN 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'lssfull0'
          par.random_sample = -1
      END 
      4: BEGIN 
          ;; NO EDGE CUTS!! 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'lssfull0'
          par.random_sample = -1
          par.max_allowed_angle = 20.0
          par.edgecuts = 0

      END 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; vagc with redshift cuts
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      5: BEGIN 
          rmax = 11500.0
          nbin = 18
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'vagc_zcut'
          par.source_sample = 6 ; r > 1/3
          par.random_sample = 9 ; The 10 Mpc sample
      END 
      6: BEGIN 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'vagc_zcut'
          par.source_sample = 6 ; r > 1/3
          par.random_sample = 9 
      END 


      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; lowz
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      7: BEGIN 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'lowz'
          par.source_sample = 6 ; r > 1/3
          par.max_allowed_angle = 20.0
          par.random_sample = 3
      END 
      8: BEGIN 
          rmax = 2801.8093
          nbin = 14
          sigmacrit_style = 1
          shape_correction_style = 3 ; princeton regauss
          self.catalog = 'lowz'
          par.source_sample = 6 ; r > 1/3
          par.random_sample = 9 

          ;; No edge cuts!!!
          par.max_allowed_angle = 20.0
          par.edgecuts = 0
      END 

      ELSE: message,'Unknown sample: '+ntostr(sample)
  ENDCASE 

  self.sample = sample
  par.catalog = self.catalog

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


END 

FUNCTION vagc_lensing::catalog
  return, self.catalog
END 




;; This has minimal info; do we just want to 
;; join on vagc later, or create a catalog with
;; more stuff in it?
FUNCTION vagc_lensing::zrange, string=string
  IF keyword_set(string) THEN BEGIN 
      return,['0.01','0.4']
  ENDIF ELSE BEGIN 
      return,[0.01, 0.4]
  ENDELSE 
END 
function vagc_lensing::get, columns=columns, status=status

    case self->catalog() of
        ;; Well-understood for lss
        'lssfull0': begin 
            query = $
                'SELECT '+$
                'l.vagc_index, '+$
                'l.ra, l.dec, l.z, l.fgot, '+$
                'v.modelabsmag as absmag, v.kmodelfratio as kfratio, '+$
                'v.vdisp '+$
                'FROM '+$
                'lssfull0 as l, vagc as v '+$
                'WHERE '+$
                'v.vagc_index = l.vagc_index'

        end 
        'vagc_zcut': begin 
            ;; everything at decent redshift
            zrange = self->zrange(/string)
            zmin = zrange[0]
            zmax = zrange[1]
            if n_elements(columns) ne 0 then begin 
                cols = strjoin(columns, ',')
            endif else begin 
                cols ='vagc_index, ra, dec, z, modelabsmag as absmag'
            endelse 
            query = $
                'SELECT '+$
                cols+' '+$
                'FROM '+$
                'vagc '+$
                'WHERE '+$
                'z > '+zmin+' AND z < '+zmax
        end 
        'lowz': begin 
            if n_elements(columns) ne 0 then begin 
                cols = strjoin(columns, ',')
            endif else begin 
                cols ='vagc_index, ra, dec, zdist as z, absmag, rc3dist '
            endelse 
            query = $
                'SELECT '+$
                    cols+' '+$
                'FROM '+$
                    'lowz'
        end 
        else: message,'Unsupported catalog: '+self->catalog()
    endcase 
 
    print
    print,query
    struct = self->postgres::query(query, status=status)
    if status ne self->postgres::status_val('success') then begin 
        message,'Query failed to return results. ',/inf
        message,'Halting'
    endif 
    return, struct
end 






; From Michael - I'd divide them into three classes (low, medium, "high"
; luminosity, far away from bright objects):
;  ilow=where(comb.absmag[2] gt -16. and comb.absmag[2] lt -13. and $
;             comb.rc3dist gt 0.5)
;  imedium=where(comb.absmag[2] gt -17.5 and comb.absmag[2] lt -16. and $
;                comb.rc3dist gt 0.5)
;  ihigh=where(comb.absmag[2] gt -19. and comb.absmag[2] lt -17.5 and $
;              comb.rc3dist gt 0.5)

;; r-band absmag
PRO vagc_lensing::absmagr_bins, nbin, lowlim, highlim

  CASE nbin OF
      3: BEGIN 
          IF self->catalog() NE 'lowz' THEN BEGIN 
              message,'Must be lowz catalog for absmagr3'
          ENDIF 
          ;; lowlim  = [-19.0,-17.5,-16.0]
          ;; highlim = [-17.5,-16.0,-13.0]

          lowlim  = [-16.0, -17.5, -19.0]
          highlim = [-13.0, -16.0, -17.5]
      END 
      ELSE: message,'Unsupported nbin: '+ntostr(nbin)
  ENDCASE 

END 

FUNCTION vagc_lensing::absmagr_where_string, nbin, labels=labels

  self->absmagr_bins, nbin, lowlim, highlim
  nstr = 'M!dr!N'
  FOR i=0L, nbin-1 DO BEGIN 

      lowstr = ntostr(lowlim[i],5,/round)
      highstr = ntostr(highlim[i],5,/round)

      tlabel = lowstr+' < '+nstr+' < '+highstr
      add_arrval, tlabel, labels

      tstring = $
        'struct.absmag[2] GE '+ntostr(lowlim[i])+' AND '+$
        'struct.absmag[2] LT '+ntostr(highlim[i])+' AND '+$
        'struct.rc3dist GT 0.5'
      add_arrval, tstring, where_string

  ENDFOR 
  return, where_string
END 

FUNCTION vagc_lensing::subtype_nbin, subtype
  IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: nbin=mb->nbin(subtype)'
  ENDIF 
  CASE subtype OF
      'absmagr3': return,3
      ELSE: message,'Unknown subtyep: '+strn(subtype)
  ENDCASE 
END 


FUNCTION vagc_lensing::where_string, subtype, nbin=nbin, labels=labels
  IF n_elements(subtype) EQ 0 THEN BEGIN 
      message,'-Syntax: ws=v->where_string(subtype)'
  ENDIF 

  nbin = self->subtype_nbin(subtype)

  cat = self->catalog()
  CASE subtype OF 
      'absmagr3': return,self->absmagr_where_string(nbin, labels=labels)
      ELSE: message,'Unknown type: '+ntostr(subtype)
  ENDCASE 

END 
FUNCTION vagc_lensing::labels, subtype, bin=bin
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





PRO vagc_lensing::combine_lensum, randnum=randnum

  CASE self->catalog() OF 
      'lowz': BEGIN 
          zbinsize = 0.0025
      END 
      ELSE: message,'Unknown catalog: '+self->catalog()
  ENDCASE 

  self->objshear::combine_lensum, randnum=randnum, zbinsize=zbinsize

END 


PRO vagc_lensing::sub, subtype, randnum=randnum


  IF n_elements(subtype) EQ 0 THEN BEGIN 
      print,'-Syntax: v->sub, subtype, randnum=, zbinsize='
      print
      message,'Halting'
  ENDIF 

  CASE self->catalog() OF 
      'lowz': BEGIN 
          average_tags = ['z', 'absmag[2]']
          zbinsize = 0.0025
      END 
      ELSE: message,'Unknown catalog: '+self->catalog()
  ENDCASE 

  ;; just to make sure we know about this subtype
  nbin = self->subtype_nbin(subtype)

  self->objshear::sub_sample, subtype, randnum=randnum, $
    average_tags=average_tags, zbinsize=zbinsize, /zhist

END 



















PRO vagc_lensing::plot_profile, type, subtype=subtype, bin=bin, $
                xlog=xlog, xmin=xmin, xmax=xmax, $
                ylog=ylog, ymin=ymin, ymax=ymax, $
                incolor=incolor, dops=dops
  
  IF n_elements(subtype) NE 0 THEN BEGIN 
      CASE subtype OF
          'absmagr3': BEGIN 
              IF n_elements(bin) EQ 0 THEN BEGIN 
                  landscape=1
                  label_charsize=0.7
                  charsize=1
              ENDIF 
              IF n_elements(xlog) EQ 0 THEN xlog=1
              IF n_elements(ylog) EQ 0 THEN ylog=1
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





PRO vagc_lensing::plot_absmag, clr, struct=struct, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      dir = self->lensdir('plot',/base)
      file = 'vagc_absmag_'+!colors[clr]+'.ps'
      file = concat_dir(dir, file)
      begplot, name=file, /landscape
  ENDIF 

  IF n_elements(struct) EQ 0 THEN BEGIN
      struct = self->getmore()
  ENDIF 

  xtitle = 'M!D'+!colors[clr]

  CASE clr OF
      0: xrange=[-8,-25]
      1: xrange=[-9,-26]
      2: xrange=[-10,-27]
      3: xrange=[-10,-27]
      4: xrange=[-10,-27]
      ELSE: message,'bad color: '+ntostr(clr)
  ENDCASE 
  plothist, struct.modelabsmag[clr], min=xrange[1], max=xrange[0], bin=0.1,$
    xrange=xrange, xstyle=3, $
    /ylog, yrange=[0.9,2.e4], ystyle=3, $
    ytickf='loglabels', $
    xtitle=xtitle, ytitle='Number'

  IF keyword_set(dops) THEN endplot, /landfix

END 





FUNCTION vagc_lensing::cleanup
  return,1
END 
PRO vagc_lensing__define
  struct = {$
             vagc_lensing, $
             sample: 0, $
             catalog: '', $
             INHERITS vagc, $   ; note vagc inherits postgres
             INHERITS objshear $
           }
END 
