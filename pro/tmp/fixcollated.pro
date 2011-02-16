PRO fixcollated, spAll, collated, newcoll, runold, rerunold, camcolold, fieldold, idold

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: fixcollated, spAll, collated, newcoll'
      return
  ENDIF 

  ;; add id tags, target tags to collated file, from spAll.fits file

  print,'Making copies of run,rerun,etc'
  key=get_kbrd(1) & IF key EQ 'q' THEN return
  runold = spAll.run
  rerunold = spAll.rerun
  camcolold = spAll.camcol
  fieldold = spAll.field
  idold = spAll.id


  ;; put run,rerun,camcol,field,id into tags
  spAll.run = spAll.objid[0]
  spAll.rerun = spAll.objid[1]
  spAll.camcol = spAll.objid[2]
  spAll.field = spAll.objid[3]
  spAll.id = spAll.objid[4]

  
  print,'Adding run,rerun,camcol,field,id,primtarget,sectarget'
  key=get_kbrd(1) & IF key EQ 'q' THEN return

  add_tags,collated, $
    ['run','rerun','camcol','field','id','primtarget','sectarget'],$
    ['0L', '0',    '0',     '0',    '0L','0L',        '0L'],$
    newcoll

  ;; put in id's
  newcoll.run = newcoll.objid[0]
  newcoll.rerun = newcoll.objid[1]
  newcoll.camcol = newcoll.objid[2]
  newcoll.field = newcoll.objid[3]
  newcoll.id = newcoll.objid[4]

  ;; remove duplicates from each, perferrably keeping the one with
  ;; primtarget set
  print,'Removing duplicates'
  key=get_kbrd(1) & IF key EQ 'q' THEN return

  rmdup_spectra, spAll, goodAll
  rmdup_spectra, newcoll, goodcoll

  ;; match
  
  print,'Matching'
  key=get_kbrd(1) & IF key EQ 'q' THEN return
  photo_match, spAll[goodAll].run, spAll[goodAll].rerun,spAll[goodAll].camcol,spAll[goodAll].field,$
               spAll[goodAll].id, $
               newcoll[goodcoll].run, newcoll[goodcoll].rerun,newcoll[goodcoll].camcol,$
               newcoll[goodcoll].field, newcoll[goodcoll].id, $
               mAll, mcoll
  mAll = goodAll[mAll]
  mcoll = goodcoll[mcoll]
  help,mAll,mcoll,newcoll

  print,'Copying primtarget'
  key=get_kbrd(1) & IF key EQ 'q' THEN return
  newcoll = temporary(newcoll[mcoll])
  newcoll.primtarget = spAll[mAll].primtarget
  newcoll.sectarget = spAll[mAll].sectarget

  spAll.run = runold
  spAll.rerun = rerunold
  spAll.camcol = camcolold
  spAll.field = fieldold
  spAll.id = idold

END 
