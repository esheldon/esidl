pro tvicl, ra, dec

    ; Just copied from Blanton's stuff, doesn't work yet need to fix
    cldir=getenv('DATA')+'/masked_clusters/'

    ihr=long(ra/15.)
    idec=long(abs(dec)/2.)*2.
    dsign='p'
    if(dec lt 0.) then dsign='m'
    outdir=getenv('DATA')+'/masked_clusters/'+string(ihr,f='(i2.2)')+'h'
    outdir=outdir+'/'+dsign+strtrim(string(idec, f='(i2.2)'),2)
    prefix='mBCG-'+hogg_iau_name(ra,dec,'')
    outdir=outdir+'/'+prefix
    rim=mrdfits(outdir[0]+'/cen-'+prefix[0]+'-r.fits.gz',0)
    gim=mrdfits(outdir[0]+'/cen-'+prefix[0]+'-g.fits.gz',0)
    iim=mrdfits(outdir[0]+'/cen-'+prefix[0]+'-i.fits.gz',0)
 
    rgbview, iim, rim, gim

    return

end
