pro ess_dbbuild,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,$
                v14,v15,v16,v17,v18,$
                 v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,$
                 v30,v31,v32,v33,v34,v35,v36,v37,v38,v39, $
                 v40,v41,v42,v43,v44,v45,v46,v47,v48,v49, $                 
                 v50,v51,v52,v53,v54,v55,v56,v57,v58,v59, $
                 v60,v61,v62,v63,v64,v65,v66,v67,v68,v69, $
                 v70,v71,v72,v73,v74,v75,v76,v77,v78,v79, $
                 v80,v81,v82,v83,v84,v85,v86,v87,v88,v89, $
                 v90,v91,v92,v93,v94,v95,v96,v97,v98,v99, $
                 v100,v101,v102,v103,v104,v105,v106,v107,v108,v109, $
                 v110,v111,v112,v113,v114,v115,v116,v117,v118,v119, $
                 v120,v121,v122,v123,v124,v125,v126,v127,v128,v129, v130, $
                 NOINDEX = noindex, STATUS=STATUS, SILENT=SILENT
;+
; NAME:
;	DBBUILD
; PURPOSE:
;	Build a database by appending new values for every item.  
; EXPLANATION:
;	The database must be opened for update (with DBOPEN) before calling 
;	DBBUILD.
;
; CALLING SEQUENCE:
;	DBBUILD, [ v1, v2, v3, v4......v130, /NOINDEX, /SILENT, STATUS =  ]
;
; INPUTS:
;	v1,v2....v130 - vectors containing values for all items in the database.
;         V1 contains values for the first item, V2 for the second, etc.
;         The number of vectors supplied must equal the number of items
;         (excluding entry number) in the database.  The number of elements 
;         in each vector should be the same.   A multiple valued item
;         should be dimensioned NVALUE by NENTRY, where NVALUE is the number
;         of values, and NENTRY is the number of entries.
;
; OPTIONAL INPUT KEYWORDS:
;	NOINDEX - If this keyword is supplied and non-zero then DBBUILD will
;             *not* create an indexed file.    Useful to save time if
;             DBBUILD is to be called several times and the indexed file need
;             only be created on the last call
;
;	SILENT  - If the keyword SILENT is set and non-zero, then DBBUILD
;	      will not print a message when the index files are generated
;
; OPTIONAL OUTPUT KEYWORD:
;	STATUS - Returns a status code denoting whether the operation was
;	      successful (1) or unsuccessful (0).  Useful when DBBUILD is
;	      called from within other applications.
;
; EXAMPLE:
;	Suppose a database named STARS contains the four items NAME,RA,DEC, and 
;	FLUX.   Assume that one already has the four vectors containing the
;	values, and that the database definition (.DBD) file already exists.
;
;	IDL> !PRIV=2                  ;Writing to database requires !PRIV=2
;	IDL> dbcreate,'stars',1,1   ;Create database (.DBF) & index (.DBX) file
;	IDL> dbopen,'stars',1         ;Open database for update
;	IDL> dbbuild,name,ra,dec,flux ;Write 4 vectors into the database
;
; NOTES:
;	Do not call DBCREATE before DBBUILD if you want to append entries to
;	an existing database
;
;	DBBUILD checks that each value vector matches the idl type given in the
;	database definition (.DBD) file, and that character strings are the 
;	proper length. 
; REVISION HISTORY:
;	Written          W. Landsman           March, 1989
;	Added /NOINDEX keyword           W. Landsman        November, 1992
;	User no longer need supply all items   W. Landsman  December, 1992 
;	Added STATUS keyword, William Thompson, GSFC, 1 April 1994
;	Added /SILENT keyword, William Thompson, GSFC, October 1995
;	Allow up to 130 items, fix problem if first item was multiple value
;				  W. Landsman    GSFC, July 1996
;	Faster build of external databases on big endian machines 
;				  W. Landsman    GSFC, November 1997  
;	Converted to IDL V5.0   W. Landsman 24-Nov-1997
;       Use SIZE(/TNAME) for error mesage display  W.Landsman   July 2001
;       Fix message display error introduced July 2001  W. Landsman   Oct. 2001 
;-
  COMPILE_OPT IDL2
  On_error,2                    ;Return to caller
  npar = N_params()
  if npar LT 1 then begin
      print,'Syntax - DBBUILD, v1, [ v2, v3, v4, v5, ... v130,' 
      print,'         /NOINDEX, /SILENT, STATUS =  ]'
      return
  endif
  
  dtype = ['UNDEFINED','BYTE','INTEGER*2','INTEGER*4','REAL*4','REAL*8', $
           'COMPLEX','STRING','STRUCTURE','DCOMPLEX','POINTER','OBJECT', $ 
           'UNSIGNED*2', 'UNSIGNED*4', 'INTEGER*8','UNSIGNED*8']
 
 
;  Initialize STATUS as unsuccessful (0).  If the routine is successful, this
;  will be updated below.
 
  status = 0
  
  nitem = db_info( 'ITEMS' )
  if N_elements( nitem ) EQ 0 then return
  
  items = indgen(nitem)
  db_item, items, itnum, ivalnum, idltype, sbyte, numvals, nbyte
  for i = 1,npar do begin     ;Get the dimensions and type of each input vector
      
      ii = strtrim(i,2)
      test = execute('s=size(v' + ii +')' )
      
      if s[s[0] + 1] NE idltype[i] then begin
          message, 'Item ' + strtrim( db_item_info('NAME',i),2) + $
            ' - parameter '+strtrim(i,2) + ' - has an incorrect data type',/INF
          message, 'Required data type is ' + dtype[idltype[i]], /INF
          message, 'Supplied data type is ' + dtype[s[s[0]+1]], /INF
          return
      endif
      
  endfor
  external = db_info('external',0)
  if external then noconvert = is_ieee_big() else noconvert = 1b
 
  maxpar = 130
  nitems = ( (npar<nitem) GE indgen(maxpar+1L))
  entry = make_array( DIMEN = db_info('LENGTH'),/BYTE ) ;Empty entry array
  nvalues = long( db_item_info( 'NVALUES' ) )       ;# of values per item
  nbyte = nbyte*nvalues                             ;Number of bytes per item
  Nv = N_elements(v1)/nvalues[1]                   
  for i = 0l, Nv - 1 do begin

       i1 = i*nvalues         
       i2 = i1 + nvalues -1
        dbxput,0l,entry,idltype[0],sbyte[0],nbyte[0]
        dbxput,v1[i1[1]:i2[1]],entry,idltype[1],sbyte[1],nbyte[1]
       if nitems[2] then begin
        dbxput,v2[i1[2]:i2[2]],entry,idltype[2],sbyte[2],nbyte[2]
       if nitems[3] then begin 
        dbxput,v3[i1[3]:i2[3]],entry,idltype[3],sbyte[3],nbyte[3]
       if nitems[4] then begin 
        dbxput,v4[i1[4]:i2[4]],entry,idltype[4],sbyte[4],nbyte[4]
       if nitems[5] then begin 
        dbxput,v5[i1[5]:i2[5]],entry,idltype[5],sbyte[5],nbyte[5]
       if nitems[6] then begin 
        dbxput,v6[i1[6]:i2[6]],entry,idltype[6],sbyte[6],nbyte[6]
       if nitems[7] then begin 
        dbxput,v7[i1[7]:i2[7]],entry,idltype[7],sbyte[7],nbyte[7]
       if nitems[8] then begin 
        dbxput,v8[i1[8]:i2[8]],entry,idltype[8],sbyte[8],nbyte[8]
       if nitems[9] then begin 
        dbxput,v9[i1[9]:i2[9]],entry,idltype[9],sbyte[9],nbyte[9]
       if nitems[10] then begin 
        dbxput,v10[i1[10]:i2[10]],entry,idltype[10],sbyte[10],nbyte[10]
       if nitems[11] then begin 
        dbxput,v11[i1[11]:i2[11]],entry,idltype[11],sbyte[11],nbyte[11]
       if nitems[12] then begin 
        dbxput,v12[i1[12]:i2[12]],entry,idltype[12],sbyte[12],nbyte[12]
       if nitems[13] then begin 
        dbxput,v13[i1[13]:i2[13]],entry,idltype[13],sbyte[13],nbyte[13]
       if nitems[14] then begin
        dbxput,v14[i1[14]:i2[14]],entry,idltype[14],sbyte[14],nbyte[14]
       if nitems[15] then begin
        dbxput,v15[i1[15]:i2[15]],entry,idltype[15],sbyte[15],nbyte[15]
       if nitems[16] then begin
        dbxput,v16[i1[16]:i2[16]],entry,idltype[16],sbyte[16],nbyte[16]
       if nitems[17] then begin
        dbxput,v17[i1[17]:i2[17]],entry,idltype[17],sbyte[17],nbyte[17]
       if nitems[18] then begin
        dbxput,v18[i1[18]:i2[18]],entry,idltype[18],sbyte[18],nbyte[18]
       if nitems[19] then begin
          dbxput,v19[i1[19]:i2[19]],entry,idltype[19],sbyte[19], nbyte[19]
       if nitems[20] then begin
          dbxput, v20[ i1[20]:i2[20] ],entry,idltype[20],sbyte[20],nbyte[20]
       if nitems[21] then begin
          dbxput,v21[i1[21]:i2[21]],entry,idltype[21],sbyte[21],nbyte[21]
       if nitems[22] then begin
          dbxput,v22[i1[22]:i2[22]],entry,idltype[22],sbyte[22],nbyte[22]
       if nitems[23] then begin
          dbxput,v23[i1[23]:i2[23]],entry,idltype[23],sbyte[23],nbyte[23]
       if nitems[24] then begin
          dbxput,v24[i1[24]:i2[24]],entry,idltype[24],sbyte[24],nbyte[24]
       if nitems[25] then $
          dbxput,v25[i1[25]:i2[25]],entry,idltype[25],sbyte[25],nbyte[25]
       if nitems[26] then $
          dbxput,v26[i1[26]:i2[26]],entry,idltype[26],sbyte[26],nbyte[26]
       if nitems[27] then $
          dbxput,v27[i1[27]:i2[27]],entry,idltype[27],sbyte[27],nbyte[27]
       if nitems[28] then $
          dbxput,v28[i1[28]:i2[28]],entry,idltype[28],sbyte[28],nbyte[28]
       if nitems[29] then $
          dbxput,v29[i1[29]:i2[29]],entry,idltype[29],sbyte[29],nbyte[29]
       if nitems[30] then $
          dbxput,v30[i1[30]:i2[30]],entry,idltype[30],sbyte[30],nbyte[30]

       if nitems[31] then $
          dbxput,v31[i1[31]:i2[31]],entry,idltype[31],sbyte[31],nbyte[31]
       if nitems[32] then $
          dbxput,v32[i1[32]:i2[32]],entry,idltype[32],sbyte[32],nbyte[32]
       if nitems[33] then $
          dbxput,v33[i1[33]:i2[33]],entry,idltype[33],sbyte[33],nbyte[33]
       if nitems[34] then $
          dbxput,v34[i1[34]:i2[34]],entry,idltype[34],sbyte[34],nbyte[34]
       if nitems[35] then $
          dbxput,v35[i1[35]:i2[35]],entry,idltype[35],sbyte[35],nbyte[35]
       if nitems[36] then $
          dbxput,v36[i1[36]:i2[36]],entry,idltype[36],sbyte[36],nbyte[36]
       if nitems[37] then $
          dbxput,v37[i1[37]:i2[37]],entry,idltype[37],sbyte[37],nbyte[37]
       if nitems[38] then $
          dbxput,v38[i1[38]:i2[38]],entry,idltype[38],sbyte[38],nbyte[38]
       if nitems[39] then $
          dbxput,v39[i1[39]:i2[39]],entry,idltype[39],sbyte[39],nbyte[39]
       if nitems[40] then $
          dbxput,v40[i1[40]:i2[40]],entry,idltype[40],sbyte[40],nbyte[40]
       if nitems[41] then $
          dbxput,v41[i1[41]:i2[41]],entry,idltype[41],sbyte[41],nbyte[41]
       if nitems[42] then $
          dbxput,v42[i1[42]:i2[42]],entry,idltype[42],sbyte[42],nbyte[42]
       if nitems[43] then $
          dbxput,v43[i1[43]:i2[43]],entry,idltype[43],sbyte[43],nbyte[43]
       if nitems[44] then $
          dbxput,v44[i1[44]:i2[44]],entry,idltype[44],sbyte[44],nbyte[44]
       if nitems[45] then $
          dbxput,v45[i1[45]:i2[45]],entry,idltype[45],sbyte[45],nbyte[45]
       if nitems[46] then $
          dbxput,v46[i1[46]:i2[46]],entry,idltype[46],sbyte[46],nbyte[46]
       if nitems[47] then $
          dbxput,v47[i1[47]:i2[47]],entry,idltype[47],sbyte[47],nbyte[47]
       if nitems[48] then $
          dbxput,v48[i1[48]:i2[48]],entry,idltype[48],sbyte[48],nbyte[48]
       if nitems[49] then $
          dbxput,v49[i1[49]:i2[49]],entry,idltype[49],sbyte[49],nbyte[49]
       if nitems[50] then $
          dbxput,v50[i1[50]:i2[50]],entry,idltype[50],sbyte[50],nbyte[50]
       if nitems[51] then $
          dbxput,v51[i1[51]:i2[51]],entry,idltype[51],sbyte[51],nbyte[51]
       if nitems[52] then $
          dbxput,v52[i1[52]:i2[52]],entry,idltype[52],sbyte[52],nbyte[52]
       if nitems[53] then $
          dbxput,v53[i1[53]:i2[53]],entry,idltype[53],sbyte[53],nbyte[53]
       if nitems[54] then $
          dbxput,v54[i1[54]:i2[54]],entry,idltype[54],sbyte[54],nbyte[54]
       if nitems[55] then $
          dbxput,v55[i1[55]:i2[55]],entry,idltype[55],sbyte[55],nbyte[55]
       if nitems[56] then $
          dbxput,v56[i1[56]:i2[56]],entry,idltype[56],sbyte[56],nbyte[56]
       if nitems[57] then $
          dbxput,v57[i1[57]:i2[57]],entry,idltype[57],sbyte[57],nbyte[57]
       if nitems[58] then $
          dbxput,v58[i1[58]:i2[58]],entry,idltype[58],sbyte[58],nbyte[58]
       if nitems[59] then $
          dbxput,v59[i1[59]:i2[59]],entry,idltype[59],sbyte[59],nbyte[59]
       if nitems[60] then $
          dbxput,v60[i1[60]:i2[60]],entry,idltype[60],sbyte[60],nbyte[60]
       if nitems[61] then $
          dbxput,v61[i1[61]:i2[61]],entry,idltype[61],sbyte[61],nbyte[61]
       if nitems[62] then $
          dbxput,v62[i1[62]:i2[62]],entry,idltype[62],sbyte[62],nbyte[62]
       if nitems[63] then $
          dbxput,v63[i1[63]:i2[63]],entry,idltype[63],sbyte[63],nbyte[63]
       if nitems[64] then $
          dbxput,v64[i1[64]:i2[64]],entry,idltype[64],sbyte[64],nbyte[64]
       if nitems[65] then $
          dbxput,v65[i1[65]:i2[65]],entry,idltype[65],sbyte[65],nbyte[65]
       if nitems[66] then $
          dbxput,v66[i1[66]:i2[66]],entry,idltype[66],sbyte[66],nbyte[66]
       if nitems[67] then $
          dbxput,v67[i1[67]:i2[67]],entry,idltype[67],sbyte[67],nbyte[67]
       if nitems[68] then $
          dbxput,v68[i1[68]:i2[68]],entry,idltype[68],sbyte[68],nbyte[68]
       if nitems[69] then $
          dbxput,v69[i1[69]:i2[69]],entry,idltype[69],sbyte[69],nbyte[69]
       if nitems[70] then $
          dbxput,v70[i1[70]:i2[70]],entry,idltype[70],sbyte[70],nbyte[70]
       if nitems[71] then $
          dbxput,v71[i1[71]:i2[71]],entry,idltype[71],sbyte[71],nbyte[71]
       if nitems[72] then $
          dbxput,v72[i1[72]:i2[72]],entry,idltype[72],sbyte[72],nbyte[72]
       if nitems[73] then $
          dbxput,v73[i1[73]:i2[73]],entry,idltype[73],sbyte[73],nbyte[73]
       if nitems[74] then $
          dbxput,v74[i1[74]:i2[74]],entry,idltype[74],sbyte[74],nbyte[74]
       if nitems[75] then $
          dbxput,v75[i1[75]:i2[75]],entry,idltype[75],sbyte[75],nbyte[75]
       if nitems[76] then $
          dbxput,v76[i1[76]:i2[76]],entry,idltype[76],sbyte[76],nbyte[76]
       if nitems[77] then $
          dbxput,v77[i1[77]:i2[77]],entry,idltype[77],sbyte[77],nbyte[77]
       if nitems[78] then $
          dbxput,v78[i1[78]:i2[78]],entry,idltype[78],sbyte[78],nbyte[78]
       if nitems[79] then $
          dbxput,v79[i1[79]:i2[79]],entry,idltype[79],sbyte[79],nbyte[79]
       if nitems[80] then $
          dbxput,v80[i1[80]:i2[80]],entry,idltype[80],sbyte[80],nbyte[80]
       if nitems[81] then $
          dbxput,v81[i1[81]:i2[81]],entry,idltype[81],sbyte[81],nbyte[81]
       if nitems[82] then $
          dbxput,v82[i1[82]:i2[82]],entry,idltype[82],sbyte[82],nbyte[82]
       if nitems[83] then $
          dbxput,v83[i1[83]:i2[83]],entry,idltype[83],sbyte[83],nbyte[83]
       if nitems[84] then $
          dbxput,v84[i1[84]:i2[84]],entry,idltype[84],sbyte[84],nbyte[84]
       if nitems[85] then $
          dbxput,v85[i1[85]:i2[85]],entry,idltype[85],sbyte[85],nbyte[85]
       if nitems[86] then $
          dbxput,v86[i1[86]:i2[86]],entry,idltype[86],sbyte[86],nbyte[86]
       if nitems[87] then $
          dbxput,v87[i1[87]:i2[87]],entry,idltype[87],sbyte[87],nbyte[87]
       if nitems[88] then $
          dbxput,v88[i1[88]:i2[88]],entry,idltype[88],sbyte[88],nbyte[88]
       if nitems[89] then $
          dbxput,v89[i1[89]:i2[89]],entry,idltype[89],sbyte[89],nbyte[89]
       if nitems[90] then $
          dbxput,v90[i1[90]:i2[90]],entry,idltype[90],sbyte[90],nbyte[90]
       if nitems[91] then $
          dbxput,v91[i1[91]:i2[91]],entry,idltype[91],sbyte[91],nbyte[91]
       if nitems[92] then $
          dbxput,v92[i1[92]:i2[92]],entry,idltype[92],sbyte[92],nbyte[92]
       if nitems[93] then $
          dbxput,v93[i1[93]:i2[93]],entry,idltype[93],sbyte[93],nbyte[93]
       if nitems[94] then $
          dbxput,v94[i1[94]:i2[94]],entry,idltype[94],sbyte[94],nbyte[94]
       if nitems[95] then $
          dbxput,v95[i1[95]:i2[95]],entry,idltype[95],sbyte[95],nbyte[95]
       if nitems[96] then $
          dbxput,v96[i1[96]:i2[96]],entry,idltype[96],sbyte[96],nbyte[96]
       if nitems[97] then $
          dbxput,v97[i1[97]:i2[97]],entry,idltype[97],sbyte[97],nbyte[97]
       if nitems[98] then $
          dbxput,v98[i1[98]:i2[98]],entry,idltype[98],sbyte[98],nbyte[98]
       if nitems[99] then $
          dbxput,v99[i1[99]:i2[99]],entry,idltype[99],sbyte[99],nbyte[99]
       if nitems[100] then $
          dbxput,v100[i1[100]:i2[100]],entry,idltype[100],sbyte[100],nbyte[100]
       if nitems[101] then $
          dbxput,v101[i1[101]:i2[101]],entry,idltype[101],sbyte[101],nbyte[101]
       if nitems[102] then $
          dbxput,v102[i1[102]:i2[102]],entry,idltype[102],sbyte[102],nbyte[102]
       if nitems[103] then $
          dbxput,v103[i1[103]:i2[103]],entry,idltype[103],sbyte[103],nbyte[103]
       if nitems[104] then $
          dbxput,v104[i1[104]:i2[104]],entry,idltype[104],sbyte[104],nbyte[104]
       if nitems[105] then $
          dbxput,v105[i1[105]:i2[105]],entry,idltype[105],sbyte[105],nbyte[105]
       if nitems[106] then $
          dbxput,v106[i1[106]:i2[106]],entry,idltype[106],sbyte[106],nbyte[106]
       if nitems[107] then $
          dbxput,v107[i1[107]:i2[107]],entry,idltype[107],sbyte[107],nbyte[107]
       if nitems[108] then $
          dbxput,v108[i1[108]:i2[108]],entry,idltype[108],sbyte[108],nbyte[108]
       if nitems[109] then $
          dbxput,v109[i1[109]:i2[109]],entry,idltype[109],sbyte[109],nbyte[109]
       if nitems[110] then $
          dbxput,v110[i1[110]:i2[110]],entry,idltype[110],sbyte[110],nbyte[110]
       if nitems[111] then $
          dbxput,v111[i1[111]:i2[111]],entry,idltype[111],sbyte[111],nbyte[111]
       if nitems[112] then $
          dbxput,v112[i1[112]:i2[112]],entry,idltype[112],sbyte[112],nbyte[112]
       if nitems[113] then $
          dbxput,v113[i1[113]:i2[113]],entry,idltype[113],sbyte[113],nbyte[113]
       if nitems[114] then $
          dbxput,v114[i1[114]:i2[114]],entry,idltype[114],sbyte[114],nbyte[114]
       if nitems[115] then $
          dbxput,v115[i1[115]:i2[115]],entry,idltype[115],sbyte[115],nbyte[115]
       if nitems[116] then $
          dbxput,v116[i1[116]:i2[116]],entry,idltype[116],sbyte[116],nbyte[116]
       if nitems[117] then $
          dbxput,v117[i1[117]:i2[117]],entry,idltype[117],sbyte[117],nbyte[117]
       if nitems[118] then $
          dbxput,v118[i1[118]:i2[118]],entry,idltype[118],sbyte[118],nbyte[118]
       if nitems[119] then $
          dbxput,v119[i1[119]:i2[119]],entry,idltype[119],sbyte[119],nbyte[119]
       if nitems[120] then $
          dbxput,v120[i1[120]:i2[120]],entry,idltype[120],sbyte[120],nbyte[120]
       if nitems[121] then $
          dbxput,v121[i1[121]:i2[121]],entry,idltype[121],sbyte[121],nbyte[121]
       if nitems[122] then $
          dbxput,v122[i1[122]:i2[122]],entry,idltype[122],sbyte[122],nbyte[122]
       if nitems[123] then $
          dbxput,v123[i1[123]:i2[123]],entry,idltype[123],sbyte[123],nbyte[123]
       if nitems[124] then $
          dbxput,v124[i1[124]:i2[124]],entry,idltype[124],sbyte[124],nbyte[124]
       if nitems[125] then $
          dbxput,v125[i1[125]:i2[125]],entry,idltype[125],sbyte[125],nbyte[125]
       if nitems[126] then $
          dbxput,v126[i1[126]:i2[126]],entry,idltype[126],sbyte[126],nbyte[126]
       if nitems[127] then $
          dbxput,v127[i1[127]:i2[127]],entry,idltype[127],sbyte[127],nbyte[127]
       if nitems[128] then $
          dbxput,v128[i1[128]:i2[128]],entry,idltype[128],sbyte[128],nbyte[128]
       if nitems[129] then $
          dbxput,v129[i1[129]:i2[129]],entry,idltype[129],sbyte[129],nbyte[129]
       if nitems[130] then $
          dbxput,v130[i1[130]:i2[130]],entry,idltype[130],sbyte[130],nbyte[130]


     endif & endif & endif & endif & endif & endif & endif & endif & endif
     endif & endif & endif & endif & endif & endif & endif & endif & endif
     endif & endif & endif & endif & endif

     dbwrt,entry,noconvert=noconvert        ;Write the entry into the database

  endfor

  if not keyword_set( NOINDEX ) then begin

      indexed = db_item_info( 'INDEX' )      ;Need to create an indexed file?
      if total(indexed) GE 1 then begin
	   if not keyword_set(silent) then	$
	           message,'Now creating indexed files',/INF
           dbindex,items
       endif

  endif

  dbclose

;  Mark successful completion, and return.

  status = 1
  return
  end
