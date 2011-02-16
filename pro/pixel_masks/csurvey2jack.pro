FUNCTION csurvey2jack, clambda, ceta, file=file, jack=jack, pixelnums=pixelnums

	IF n_params() LT 2 THEN BEGIN 
		print,'-Syntax: jackknife_id = csurvey2jack(clambda, ceta, pixelnums=, jack=)'
		return,-1
	ENDIF 

	nlam = n_elements(clambda) & neta = n_elements(ceta)
	IF nlam NE neta THEN message,'clambda and ceta must be same size'

	IF n_elements(jack) EQ 0 THEN BEGIN 

		;; This is sorted by pixelnum
		IF n_elements(file) EQ 0 THEN BEGIN 
			file = esheldon_config('jackknife_file')
		ENDIF 

		print
		print,'Reading jackknife file: ',file
		jack = read_idlstruct(file,status=status)
		IF status NE 0 THEN message,'Failed to read file'
	ENDIF 

	print
	print,'Converting clambda/ceta to pixelnum'
	pixelnums = csurvey2pix(clambda, ceta, resolution=256)

	print
	print,'Matching jackknife samples to pixel numbers'

	jackknife_ids = replicate(-1, nlam)

	minjackpix = min(jack.pixelnum, max=maxjackpix)
	mininputpix = min(pixelnums, max=maxinputpix)
	minpix = min( [minjackpix, mininputpix] )
	maxpix = max( [maxjackpix, maxinputpix] )

	hjack  = histogram(jack.pixelnum - minpix, $
		min=0, $
		max = maxpix-minpix, rev = revjack)
	hinput = histogram(pixelnums - minpix, $
		min=0, $
		max=maxpix-minpix, rev = revinput)

	njack = n_elements(jack)

	FOR i=0L, njack-1 DO BEGIN 

		pix = jack[i].pixelnum - minpix

		IF revjack[pix] NE revjack[pix+1] THEN BEGIN 


			IF revinput[pix] NE revinput[pix+1] THEN BEGIN 

				; wjack should be a single element
				wjack = revjack[ revjack[pix]:revjack[pix+1] -1 ]

				w = revinput[ revinput[pix]:revinput[pix+1] -1 ]
				jackknife_ids[w] = jack[wjack].jackknife_sample

			ENDIF 


		ENDIF 

	ENDFOR 
	return,jackknife_ids

END 
