PRO shapelet_fix_structs

  ;; Fix the outputs, which didn't include enough info

  stripes = [9,10,11,12,13,14,15,16, $
             26,27,28,29,30,31,32,33,34,35,36,37, $
             42,43,44,$
             76,82,86]

  nstripe = n_elements(stripes) 
  indir1 = '~/shapelet_outputs/'
  indir2 = vagc_lensinput_dir()
  FOR i=0L, nstripe-1 DO BEGIN 
      
      sstr = stripe2string(stripes[i])

      infile1 = indir1 + 'stripe'+sstr+'_vagc_shapelets_nmax15.fit'
      infile2 = indir2 + 'stripe'+sstr+'_vagc_matchlocal.fit'

      outFile = indir1 + 'stripe'+sstr+'-vagc-shapelets-approx-nmax15.fit'

      print
      print,infile1
      shst = mrdfits(infile1,1)

      print,infile2
      mst = mrdfits(infile2,1)
      mst.rerun = mst.match_rerun
      mst.id = mst.match_id
      

      addstruct = { $
                    z:0.0, $
                    rowc:fltarr(5), colc:fltarr(5), $
                    m_e1:fltarr(5), m_e2:fltarr(5), $
                    m_e1_psf:fltarr(5), m_e2_psf:fltarr(5), $
                    m_rr_cc:fltarr(5), m_rr_cc_psf:fltarr(5) $
                  }
      outst = create_struct(shst[0], addstruct)

      outst = replicate(outst, n_elements(shst))

      copy_struct, shst, outst

      outst.ra = mst[shst.vagc_index].ra
      outst.dec = mst[shst.vagc_index].dec

      outst.z = mst[shst.vagc_index].z

      outst.rowc = mst[shst.vagc_index].rowc
      outst.colc = mst[shst.vagc_index].colc

      outst.m_rr_cc = mst[shst.vagc_index].m_rr_cc
      outst.m_rr_cc_psf = mst[shst.vagc_index].m_rr_cc_psf

      outst.m_e1 = mst[shst.vagc_index].m_e1
      outst.m_e2 = mst[shst.vagc_index].m_e2

      outst.m_e1_psf = mst[shst.vagc_index].m_e1_psf
      outst.m_e2_psf = mst[shst.vagc_index].m_e2_psf

      ;; Fix old defaults
      w=where(shst.scale EQ 0, nw)
      IF nw NE 0 THEN BEGIN 
          print
          print,'FIXING OLD DEFAULTS'
          shst[w].scale = -1
          shst[w].psf_scale = -1
          shst[w].ortho = -1
      ENDIF 

      print,'Writing output file: ',outfile
      mwrfits, outst, outFile, /create

  ENDFOR 
END 
