pro make_admomtags,taglist, default=default, addphotomom=addphotomom

  if n_params() eq 0 then begin
    print,'-syntax make_admomtags, taglist, default=default'
    return
  endif
  
  taglist=['ID',              $
           'PARENT',          $
           'NCHILD',          $
           'OBJC_TYPE',       $
           'TYPE',            $
           'FLAGS',           $
           'FLAGS2',          $
           'OBJC_FLAGS',      $
           'OBJC_FLAGS2',     $
           'OBJC_ROWC',       $
           'OBJC_COLC',       $
           'ROWC',            $
           'COLC',            $
           'COUNTS_MODEL',    $
           'COUNTS_MODELERR', $
           'COUNTS_EXP',      $
           'COUNTS_EXPERR',   $
           'PETROCOUNTS',     $
           'PETROCOUNTSERR',  $
           'PETRORAD',        $
           'PETRORADERR',     $
           'PETROR50',        $
           'PETROR50ERR',     $
           'PETROR90',        $
           'PETROR90ERR',     $
           'PSFCOUNTS',       $
           'PSFCOUNTSERR',    $
           'STATUS',          $
           'RA',              $
           'DEC',             $
           'PRIMTARGET',      $
           'SECTARGET',       $
           'REDDENING'        ]

  IF keyword_set(addphotomom) THEN BEGIN 
      taglist = [taglist, $
                 'm_rr_cc',     $
                 'm_rr_ccerr',  $
                 'm_cr4',       $
                 'm_rr_cc_psf', $
                 'm_cr4_psf',   $
                 'm_e1',        $
                 'm_e2',        $
                 'm_e1e1err',   $
                 'm_e1e2err',   $
                 'm_e2e2err',   $
                 'm_e1_psf',    $
                 'm_e2_psf'     ]
  ENDIF 

  IF NOT keyword_set(default) THEN BEGIN
      addtags = [$
                  'SEEING',         $
                  'SKY',            $
                  'SKYERR',         $
                  'STARFLAG',       $
                  'IXX',            $
                  'IYY',            $
                  'IXY',            $
                  'RHO4',           $
                  'WHYFLAG',        $
                  'PSFIXX',         $
                  'PSFIYY',         $
                  'PSFIXY',         $
                  'PSFRHO4',        $
                  'E1',             $
                  'E2',             $
                  'MOMERR',         $
                  'ROTATION',       $
                  'R']
      taglist = [taglist, addtags]
  ENDIF

  return
end
