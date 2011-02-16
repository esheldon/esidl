PRO read_sxtsgals, tsgals
  
  file='/sdss2/data0/sdss/sxoutput/tsgals_fnalchunk.dat'
  outf='/sdss2/data0/sdss/sxoutput/tsgals_fnalchunk.fit'

  ;; if we already  made the fits file just read it in
  IF fexist(outf) THEN BEGIN
      tsgals=mrdfits(outf,1)
      return
  ENDIF 

  arrval=fltarr(5)

  tmpstruct1=create_struct('run', 0L, $
                           'rerun', 0, $
                           'camcol', 0, $
                           'field', 0, $
                           'id', 0L, $
                           'objc_rowc', 0.0,$
                           'objc_colc', 0.0, $
                           'rowc0',0.,'rowc1',0.,'rowc2',0.,'rowc3',0.,'rowc4',0.,$
                           'colc0',0.,'colc1',0.,'colc2',0.,'colc3',0.,'colc4',0.)

  tmpstruct2=create_struct('cm0',0.,'cm1',0.,'cm2',0.,'cm3',0.,'cm4',0.,$
                           'cme0',0.,'cme1',0.,'cme2',0.,'cme3',0.,'cme4',0.,$
                           'petc0',0.,'petc1',0.,'petc2',0.,'petc3',0.,'petc4',0.)
  
  tmpstruct3=create_struct('petce0',0.,'petce1',0.,'petce2',0.,'petce3',0.,'petce4',0.,$
                           'petr0',0.,'petr1',0.,'petr2',0.,'petr3',0.,'petr4',0.,$
                           'petre0',0.,'petre1',0.,'petre2',0.,'petre3',0.,'petre4',0.,$
                           'petr500',0.,'petr501',0.,'petr502',0.,'petr503',0.,'petr504',0.)

  tmpstruct4=create_struct('petr50e0',0.,'petr50e1',0.,'petr50e2',0.,'petr50e3',0.,'petr50e4',0.,$
                           'petr900',0.,'petr901',0.,'petr902',0.,'petr903',0.,'petr904',0.,$
                           'petr90e0',0.,'petr90e1',0.,'petr90e2',0.,'petr90e3',0.,'petr90e4',0.,$
                           'primtarget',0L,$
                           'sectarget',0L,$
                           'red0',0.,'red1',0.,'red2',0.,'red3',0.,'red4',0.,$
                           'ra',0d,$
                           'dec',0d)

  tmpstruct=create_struct(tmpstruct1,tmpstruct2)
  tmpstruct=create_struct(tmpstruct,tmpstruct3)
  tmpstruct=create_struct(tmpstruct,tmpstruct4)

  arrval=fltarr(5)
  tsstruct=create_struct('run', 0L, $
                         'rerun', 0, $
                         'camcol', 0, $
                         'field', 0, $
                         'id', 0L, $
                         'objc_rowc', 0.0,$
                         'objc_colc', 0.0, $
                         'rowc', arrval,$
                         'colc', arrval,$
                         'counts_model', arrval,$
                         'counts_modelerr',arrval,$
                         'petrocounts',arrval,$
                         'petrocountserr',arrval,$
                         'petrorad',arrval,$
                         'petroraderr',arrval,$
                         'petror50',arrval,$
                         'petror50err',arrval,$
                         'petror90',arrval,$
                         'petror90err',arrval,$
                         'primtarget',0L,$
                         'sectarget',0L,$
                         'reddening',arrval,$
                         'ra',0d,$
                         'dec',0d)


  rdfloatstr, file, tmpstruct, outstruct, /double

  nobj=n_elements(outstruct)

  tsgals=replicate(tsstruct,nobj)

  tsgals.run=outstruct.run
  tsgals.rerun=outstruct.rerun
  tsgals.camcol=outstruct.camcol
  tsgals.field=outstruct.field
  tsgals.id=outstruct.id
  tsgals.objc_rowc=outstruct.objc_rowc
  tsgals.objc_colc=outstruct.objc_colc

  tsgals.rowc[0] = outstruct.rowc0
  tsgals.rowc[1] = outstruct.rowc1
  tsgals.rowc[2] = outstruct.rowc2
  tsgals.rowc[3] = outstruct.rowc3
  tsgals.rowc[4] = outstruct.rowc4

  tsgals.colc[0] = outstruct.colc0
  tsgals.colc[1] = outstruct.colc1
  tsgals.colc[2] = outstruct.colc2
  tsgals.colc[3] = outstruct.colc3
  tsgals.colc[4] = outstruct.colc4

  tsgals.counts_model[0] = outstruct.cm0
  tsgals.counts_model[1] = outstruct.cm1
  tsgals.counts_model[2] = outstruct.cm2
  tsgals.counts_model[3] = outstruct.cm3
  tsgals.counts_model[4] = outstruct.cm4

  tsgals.counts_modelerr[0] = outstruct.cme0
  tsgals.counts_modelerr[1] = outstruct.cme1
  tsgals.counts_modelerr[2] = outstruct.cme2
  tsgals.counts_modelerr[3] = outstruct.cme3
  tsgals.counts_modelerr[4] = outstruct.cme4

  tsgals.petrocounts[0] = outstruct.petc0
  tsgals.petrocounts[1] = outstruct.petc1
  tsgals.petrocounts[2] = outstruct.petc2
  tsgals.petrocounts[3] = outstruct.petc3
  tsgals.petrocounts[4] = outstruct.petc4

  tsgals.petrocountserr[0] = outstruct.petce0
  tsgals.petrocountserr[1] = outstruct.petce1
  tsgals.petrocountserr[2] = outstruct.petce2
  tsgals.petrocountserr[3] = outstruct.petce3
  tsgals.petrocountserr[4] = outstruct.petce4

  tsgals.petrorad[0] = outstruct.petr0
  tsgals.petrorad[1] = outstruct.petr1
  tsgals.petrorad[2] = outstruct.petr2
  tsgals.petrorad[3] = outstruct.petr3
  tsgals.petrorad[4] = outstruct.petr4

  tsgals.petroraderr[0] = outstruct.petre0
  tsgals.petroraderr[1] = outstruct.petre1
  tsgals.petroraderr[2] = outstruct.petre2
  tsgals.petroraderr[3] = outstruct.petre3
  tsgals.petroraderr[4] = outstruct.petre4

  tsgals.petror50[0] = outstruct.petr500
  tsgals.petror50[1] = outstruct.petr501
  tsgals.petror50[2] = outstruct.petr502
  tsgals.petror50[3] = outstruct.petr503
  tsgals.petror50[4] = outstruct.petr504

  tsgals.petror50err[0] = outstruct.petr50e0
  tsgals.petror50err[1] = outstruct.petr50e1
  tsgals.petror50err[2] = outstruct.petr50e2
  tsgals.petror50err[3] = outstruct.petr50e3
  tsgals.petror50err[4] = outstruct.petr50e4

  tsgals.petror90[0] = outstruct.petr900
  tsgals.petror90[1] = outstruct.petr901
  tsgals.petror90[2] = outstruct.petr902
  tsgals.petror90[3] = outstruct.petr903
  tsgals.petror90[4] = outstruct.petr904

  tsgals.petror90err[0] = outstruct.petr90e0
  tsgals.petror90err[1] = outstruct.petr90e1
  tsgals.petror90err[2] = outstruct.petr90e2
  tsgals.petror90err[3] = outstruct.petr90e3
  tsgals.petror90err[4] = outstruct.petr90e4

  tsgals.primtarget=outstruct.primtarget
  tsgals.sectarget=outstruct.sectarget

  tsgals.reddening[0] = outstruct.red0
  tsgals.reddening[1] = outstruct.red1
  tsgals.reddening[2] = outstruct.red2
  tsgals.reddening[3] = outstruct.red3
  tsgals.reddening[4] = outstruct.red4

  tsgals.ra=outstruct.ra
  tsgals.dec=outstruct.dec

  mwrfits, tsgals, outf, /create

END 
