
PRO get_fatsoruns, printout=printout

  ;; create this file on fatso with 
  ;; print_struct,!run_status,file='fatso_runlist.txt'

  dir = '/net/cheops1/data0/imaging/dbm/'
  fatsolist = dir+'fatso_runlist.txt'
  getfromfermi = dir + 'getfromfermi.txt'

  readcol, fatsolist, $
           fnum,fstripes,fstrips,fruns,freruns,fts,fatl,fadat,fbad,fbad2,$
           format='I,I,A,I,L,F,F,F,L,L'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find runs for which we already have tsObj
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  make_runstatus_struct, rs
  rs.tsobj_exist = 'Y'
  runstatus_select, rs, si

  runs   = !run_status[si].run
  reruns = !run_status[si].rerun
  stripes = !run_status[si].stripe
  strips = !run_status[si].strip

  ;; create indices for each
  ten = ulong64(10)
  findex = ulong64(freruns)*ten^13 + ulong64(fruns)*ten^16

  index = ulong64(reruns)*ten^13 + ulong64(runs)*ten^16

  match, index, findex, m, fm, /sort

  remove, m, runs, reruns, stripes, strips
  remove, fm, fruns, freruns, fstripes, fstrips

  ;; print those that fatso has but we
  ;; don't

  IF keyword_set(printout) THEN BEGIN 
      print,'Outputting to file: '+getfromfermi
      openw, lun, getfromfermi, /get_lun
  ENDIF ELSE lun=-1

  key = ntostr(fstripes)+'-'+fstrips
  s=sort(key)
  printf,lun
  printf,lun,'runs at pitt that did not match to cheops'
  colprint,key[s],fruns[s],freruns[s], lun=lun

  key = ntostr(stripes)+'-'+strips
  s=sort(key)
  printf,lun
  printf,lun,'runs at cheops that did not match to pitt'
  colprint,key[s],runs[s],reruns[s], lun=lun

  IF keyword_set(printout) THEN free_lun, lun



END 
