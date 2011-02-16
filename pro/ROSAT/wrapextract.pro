PRO wrapextract

;;; Note:  745 and up have different number of tags so read_photo_col will
;;; crash if using default structure name.

runlist = [752, 756]
n=n_elements(runlist)

FOR i=0, n-1 DO BEGIN

  rosat_extract, runlist[i]

ENDFOR

return
end
