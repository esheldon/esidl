PRO read_sxspecgals, specgals

  file='/sdss2/data0/sdss/sxoutput/specgals.dat'
  outf='/sdss2/data0/sdss/sxoutput/specgals.fit'

  ;; if we already  made the fits file just read it in
  IF fexist(outf) THEN BEGIN
      specgals=mrdfits(outf,1)
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
                          'dec',0d,$
                          'z1d',0.0,$
                          'z1derr',0.0,$
                          'zconf',0.0,$
                          'specclass',0,$
                          'veldisp',0.0,$
                          'veldisperr',0.0)

  tmpstruct=create_struct(tmpstruct1,tmpstruct2)
  tmpstruct=create_struct(tmpstruct,tmpstruct3)
  tmpstruct=create_struct(tmpstruct,tmpstruct4)


  specstruct=create_struct('run', 0L, $
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
                           'dec',0d,$
                           'z1d',0.0,$
                           'z1derr',0.0,$
                           'zconf',0.0,$
                           'specclass',0,$
                           'veldisp',0.0,$
                           'veldisperr',0.0)


  rdfloatstr, file, tmpstruct, outstruct, /double

  nobj=n_elements(outstruct)

  specgals=replicate(specstruct,nobj)

  specgals.run=outstruct.run
  specgals.rerun=outstruct.rerun
  specgals.camcol=outstruct.camcol
  specgals.field=outstruct.field
  specgals.id=outstruct.id
  specgals.objc_rowc=outstruct.objc_rowc
  specgals.objc_colc=outstruct.objc_colc

  specgals.rowc[0] = outstruct.rowc0
  specgals.rowc[1] = outstruct.rowc1
  specgals.rowc[2] = outstruct.rowc2
  specgals.rowc[3] = outstruct.rowc3
  specgals.rowc[4] = outstruct.rowc4

  specgals.colc[0] = outstruct.colc0
  specgals.colc[1] = outstruct.colc1
  specgals.colc[2] = outstruct.colc2
  specgals.colc[3] = outstruct.colc3
  specgals.colc[4] = outstruct.colc4

  specgals.counts_model[0] = outstruct.cm0
  specgals.counts_model[1] = outstruct.cm1
  specgals.counts_model[2] = outstruct.cm2
  specgals.counts_model[3] = outstruct.cm3
  specgals.counts_model[4] = outstruct.cm4

  specgals.counts_modelerr[0] = outstruct.cme0
  specgals.counts_modelerr[1] = outstruct.cme1
  specgals.counts_modelerr[2] = outstruct.cme2
  specgals.counts_modelerr[3] = outstruct.cme3
  specgals.counts_modelerr[4] = outstruct.cme4

  specgals.petrocounts[0] = outstruct.petc0
  specgals.petrocounts[1] = outstruct.petc1
  specgals.petrocounts[2] = outstruct.petc2
  specgals.petrocounts[3] = outstruct.petc3
  specgals.petrocounts[4] = outstruct.petc4

  specgals.petrocountserr[0] = outstruct.petce0
  specgals.petrocountserr[1] = outstruct.petce1
  specgals.petrocountserr[2] = outstruct.petce2
  specgals.petrocountserr[3] = outstruct.petce3
  specgals.petrocountserr[4] = outstruct.petce4

  specgals.petrorad[0] = outstruct.petr0
  specgals.petrorad[1] = outstruct.petr1
  specgals.petrorad[2] = outstruct.petr2
  specgals.petrorad[3] = outstruct.petr3
  specgals.petrorad[4] = outstruct.petr4

  specgals.petroraderr[0] = outstruct.petre0
  specgals.petroraderr[1] = outstruct.petre1
  specgals.petroraderr[2] = outstruct.petre2
  specgals.petroraderr[3] = outstruct.petre3
  specgals.petroraderr[4] = outstruct.petre4

  specgals.petror50[0] = outstruct.petr500
  specgals.petror50[1] = outstruct.petr501
  specgals.petror50[2] = outstruct.petr502
  specgals.petror50[3] = outstruct.petr503
  specgals.petror50[4] = outstruct.petr504

  specgals.petror50err[0] = outstruct.petr50e0
  specgals.petror50err[1] = outstruct.petr50e1
  specgals.petror50err[2] = outstruct.petr50e2
  specgals.petror50err[3] = outstruct.petr50e3
  specgals.petror50err[4] = outstruct.petr50e4

  specgals.petror90[0] = outstruct.petr900
  specgals.petror90[1] = outstruct.petr901
  specgals.petror90[2] = outstruct.petr902
  specgals.petror90[3] = outstruct.petr903
  specgals.petror90[4] = outstruct.petr904

  specgals.petror90err[0] = outstruct.petr90e0
  specgals.petror90err[1] = outstruct.petr90e1
  specgals.petror90err[2] = outstruct.petr90e2
  specgals.petror90err[3] = outstruct.petr90e3
  specgals.petror90err[4] = outstruct.petr90e4

  specgals.primtarget=outstruct.primtarget
  specgals.sectarget=outstruct.sectarget

  specgals.reddening[0] = outstruct.red0
  specgals.reddening[1] = outstruct.red1
  specgals.reddening[2] = outstruct.red2
  specgals.reddening[3] = outstruct.red3
  specgals.reddening[4] = outstruct.red4

  specgals.ra=outstruct.ra
  specgals.dec=outstruct.dec
  specgals.z1d=outstruct.z1d
  specgals.z1derr=outstruct.z1derr
  specgals.zconf=outstruct.zconf
  specgals.specClass=outstruct.specClass
  specgals.velDisp=outstruct.velDisp
  specgals.velDisperr=outstruct.velDisperr

  mwrfits, specgals, outf, /create

END 
