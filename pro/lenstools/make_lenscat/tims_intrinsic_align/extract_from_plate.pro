pro extract_from_plate,pl,pln,keepnomatch=keepnomatch

    rem_tags=['CATID','OBJC_ROWC',$
	'OBJC_ROWCERR','OBJC_COLC','OBJC_COLCERR','ROWV','ROWVERR',$
	'COLV','COLVERR','SKY','SKYERR',$
	'FIBERCOUNTS','FIBERCOUNTSERR','Q','QERR','U','UERR','ISO_ROWC',$
	'ISO_ROWCERR','ISO_ROWCGRAD','ISO_COLC','ISO_COLCERR','ISO_COLCGRAD',$
	'ISO_A','ISO_AERR','ISO_AGRAD','ISO_B','ISO_BERR','ISO_BGRAD',$
	'ISO_PHI','ISO_PHIERR','ISO_PHIGRAD','R_DEV','R_DEVERR','AB_DEV',$
	'AB_DEVERR','PHI_DEV','PHI_DEVERR','COUNTS_DEV','COUNTS_DEVERR',$
	'R_EXP','R_EXPERR','AB_EXP','AB_EXPERR','PHI_EXP','PHI_EXPERR',$
	'COUNTS_EXP','COUNTS_EXPERR','TEXTURE','STAR_L','EXP_L','DEV_L',$
	'FRACPSF','TYPE','NPROF','PROFMEAN','PROFERR',$
	'STATUS','OFFSETRA','OFFSETDEC','PROPERMOTIONMATCH',$
	'PROPERMOTIONDELTA','PROPERMOTION','PROPERMOTIONANGLE','USNOBLUE',$
	'USNORED','FIRSTMATCH','FIRSTID','FIRSTLAMBDA','FIRSTETA',$
	'FIRSTDELTA','FIRSTPEAK','FIRSTINT','FIRSTRMS','FIRSTMAJOR',$
	'FIRSTMINOR','FIRSTPA','ROSATMATCH','ROSATDELTA','ROSATPOSERR',$
	'ROSATCPS','ROSATCPSERR','ROSATHR1','ROSATHR1ERR','ROSATHR2',$
	'ROSATHR2ERR','ROSATEXT','ROSATEXTLIKE','ROSATDETECTLIKE',$
	'ROSATEXPOSURE','PRIORITY','MATCHID','PLATE','FIBERID','OBJID',$
	'DIAG_PRIMTARGET','DIAG_SECTARGET','XFOCAL','YFOCAL','DIAG_RA',$
	'DIAG_DEC','SKYRES','UGRIZFIBRE','GRIFLUX','GRISN','DERIV2','MEAN',$
	'CONTCHI2','BLUESLOPE','REDSLOPE','CLASS','PLUGMAPOBJ_OBJID',$
	'HOLETYPE','PLUGMAPOBJ_RA','PLUGMAPOBJ_DEC','MAG','STARL','EXPL',$
	'DEVAUCL','OBJTYPE','PLUGMAPOBJ_XFOCAL','PLUGMAPOBJ_YFOCAL',$
	'SPECTROGRAPHID','PLUGMAPOBJ_FIBERID','THROUGHPUT',$
	'PLUGMAPOBJ_PRIMTARGET','PLUGMAPOBJ_SECTARGET','USECLOSEST']

    pln = remove_tags(pl,rem_tags)

      ;; select objects targeted as galaxies and matched to
      ;; adatc files, plus some additional cuts

      edef = 0.0
      IF keyword_set(keepnomatch) THEN BEGIN 
          j=where(pln.z gt 0.01 and pln.z lt 0.5 and $
                  pln.zwarning eq 0 and $
                  pln.petrocounts(2) gt 12 and pln.petrocounts(2) lt 20 and $
                  (pln.primtarget eq 32 or pln.primtarget eq 64 or $
                   pln.primtarget eq 96))
      ENDIF ELSE BEGIN 
          j=where(pln.matchrerun NE -1 AND $
                  pln.z gt 0.01 and pln.z lt 0.5 and $
                  pln.zwarning eq 0 and $
                  pln.petrocounts(2) gt 12 and pln.petrocounts(2) lt 20 and $
                  (pln.primtarget eq 32 or pln.primtarget eq 64 or $
                   pln.primtarget eq 96))
      ENDELSE 

      if (j(0) ne -1) then begin
        pln=pln(j)
      endif else begin
	pln.run(*)=0
      endelse

      return
      end
