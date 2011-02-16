;; Based on Dave's code

FUNCTION em::init
  return,1
END 



FUNCTION em::gmix, vect, ngauss, niter, guess=guess, weight=weight

  IF n_params() LT 3 THEN BEGIN 
      message,'-Syntax: gstr = em->gmix(vect, ngauss, niter, guess=, weight=)'
  ENDIF 
  
  s=size(vect)
  print,s
  if s(0) eq 2 then begin
      ndim=s(1)
      npoints=s(2)
  endif else begin
      ndim=1
      npoints=s(1)
  endelse

  if n_elements(weight) ne 0 then begin
      ww=weight
  endif else begin
      ww=replicate(1.0/npoints,npoints)
  endelse
  weight=ww

  IF n_elements(guess) NE 0then BEGIN 

      gstr = guess

  ENDIF ELSE BEGIN 

      gstr = {norm: 1.0, mean: fltarr(ndim), cov:identity(ndim)}
      gstr=replicate(gstr,ngauss)
      if ndim eq 2 then BEGIN
          mm = self->rand_init2(vect, weight)
      endif else begin
          mm = self->rand_init(ngauss, vect)
      endelse
                                ;randomly initialize the "means"
      gstr.mean = mm

  ENDELSE 



  ;; EM algorithm

  tau=fltarr(npoints,ngauss)
  beta=tau


  chi2m=replicate(100.0,npoints)
  for iter=0, niter-1 do begin

      for j=0,ngauss-1 do begin	
          like = $
            self->glike(vect, gstr[j].mean, gstr[j].cov, $
                        chi2min=chi2m,/fast)
          beta[*,j]=like
      endfor	

      if 1 then begin
          toti=beta#gstr.norm
          amat=(1.0/toti)#gstr.norm
          tau=beta*amat
      endif else begin
          for i=0,npoints-1 do begin
              tot=total(gstr.norm*beta[i,*])
              for j=0,ngauss-1 do begin
                  tau(i,j)=gstr(j).norm*beta[i,j]/tot
              endfor
          endfor
      endelse
      

      if n_elements(weight) ne 0 then begin
          ww=weight
      endif else begin
          ww=replicate(1.0/npoints,npoints)
      endelse

      gp=(tau##ww)
      gstr.norm=reform(gp)	

      www=ww##replicate(1,ndim)
      v=vect*www
      a=v#tau
      b=1.0/gstr.norm
      c=b##replicate(1,ndim)
      gstr.mean=c*a

      ;; My guess for diagonals
;      FOR j=0L, ngauss-1 DO BEGIN 

;          FOR k=0L, ndim-1 DO BEGIN 
;              gstr[j].cov[k,k] = $
;                total(tau[*,j]*ww*(vect-gstr[j].mean[k])^2)/gstr[j].norm
;          ENDFOR 
          
;      ENDFOR 


      if (0) and (iter lt niter-1) then begin
          gstr.mean(0)=gstr.mean(0)+5.0*randomu(seed,ngauss)
          gstr.mean(0)=gstr.mean(0)+5.0*randomu(seed,ngauss)
      endif
            
  endfor


  return,gstr

END 



FUNCTION em::rand_init, ngauss, vect

;give random initial positions to the gaussians
;uniformly distributed around the data range

  if n_params() LT 2 then begin
      message,'-syntax mu = em->rand_init(ngauss,vect)'
  endif

  ndim = size(vect, /n_dim)

  mu=fltarr(ndim,ngauss)
  for i=0L, ndim-1 do begin
      mn=min(vect(i,*))
      mx=max(vect(i,*))
      r=randomu(seed,ngauss)*.9+.05
      r=r*(mx-mn)+mn
      mu(i,*)=r
  endfor

  return, mu

END 

FUNCTION em::rand_init2,vect,weight
;give random initial positions to the gaussians
;uniformly distributed around the data range
;when there are only two points

  if n_params() LT 2 then begin
      message,'-syntax mu = em->rand_init2(vect, weight)'
  endif


  mu=fltarr(2,2)
  tot=total(weight)
  cenx=total(vect(0,*)*weight)/tot
  ceny=total(vect(1,*)*weight)/tot
  mu(0,0)=cenx
  mu(1,0)=ceny
  
  mm=max(weight,wmm)
  mu(0,1)=vect(0,wmm)
  mu(1,1)=vect(1,wmm)
  

  return, mu
END 






FUNCTION em::glike, x, mu, cov, prob=prob, chi2min=chi2min, fast=fast

;evaluate the likelihood at point x
;of a gaussian with mean mu and covariance cov
;x may be two dimensional like
;(ndim,npoints) and it will loop over npoints
;and then like will be of size npoints

  if n_params() LT 3 then begin
      message,'-syntax like = em->glike(x,mu,cov,prob=, chi2min=, /fast)'
  endif

  scov=size(cov)
  if scov(0) eq 0 then cov=[cov]
  invcov=invert(cov)

  s=size(x) 
  if s(0) eq 1 then BEGIN
      if n_elements(chi2min) eq 0 then chi2min=0.0
      chi2=(x-mu)##invcov##transpose(x-mu)
      chi2=chi2(0)
      chi2=chi2-chi2min
      chi2=chi2 < 40.0
      if chi2 eq 40.0 then begin
          prob=0.0 
      endif else begin
          if n_params() gt 4 then begin
              prob=1.0-chisqr_pdf(chi2,s(1))
          endif
      endelse
      like=exp(-.5*chi2)/(sqrt(determ(cov))*(2*!pi)^2)
      return,like
  endif


  if keyword_set(fast) then begin
      b=(sqrt(determ(cov))*(2*!pi)^2)
      mmu=mu#replicate(1.0,s(2))
      xx=transpose(x-mmu)
      chi2=(xx*(invcov##xx))#[1.0,1.0]
      chi2=(chi2-chi2min) < 40.0
      like=exp(-.5*chi2)/b
      return, like
  endif
  


  if s(0) eq 2 then begin
      if n_elements(chi2min) eq 0 then chi2min=replicate(0.0,s(2))
      like=fltarr(s(2))
      prob=like
      b=(sqrt(determ(cov))*(2*!pi)^2)

      for i=0, s(2)-1 do begin
          chi2=(x(*,i)-mu)##invcov##transpose(x(*,i)-mu)
          chi2=chi2(0)-chi2min(i)
          chi2=chi2 < 40.0
          if chi2 eq 40.0 then begin
              prob(i)=0.0
          endif else begin
              if n_params() gt 4 then begin
                  prob(i)=1.0-chisqr_pdf(chi2,s(1))
              endif
          endelse
          l=exp(-.5*chi2)/b
          like(i)=l
      endfor
      return,like
  endif

  print,'must be one or two dimensional'
  like=-1
  return, like
end




















PRO em::gmix_plot_model_1d, model_struct, $
  _extra=_extra, $
  overplot=overplot, nolegend=nolegend

  simpctable, colorlist=colorlist
  
  ngauss = model_struct.ngauss

  line = 2
  oline = 0
  
  IF keyword_set(overplot) THEN BEGIN 
      colors = colorlist[ 1 + lindgen(ngauss+1) ]
;      oplot, model_struct.x, model_struct.model,color=colors[0]
  ENDIF ELSE BEGIN 
      colors = colorlist[0:ngauss]
      plot, model_struct.x, model_struct.model, line=line, _extra=_extra
  ENDELSE 


  FOR i=0L, model_struct.ngauss-1 DO BEGIN

      comm = 'tmodel = model_struct.model'+ntostr(i+1)
      IF NOT execute(comm) THEN message,'Doh!'

      oplot, model_struct.x, tmodel,color=colors[i+1], line=oline
      
  ENDFOR 

  IF keyword_set(overplot) THEN BEGIN 
      oplot, model_struct.x, model_struct.model,color=colors[0], line=line
  ENDIF 

  IF NOT keyword_set(nolegend) THEN BEGIN 
      mess = ['model', 'gauss '+ntostr(1+lindgen(model_struct.ngauss))]
      legend,mess,colors=colors,box=0,line=[line, replicate(oline, ngauss)],$
             charsize=1
  ENDIF 


END 

FUNCTION em::gmix_model_1d, gstr, minx, maxx, nx

  ngauss = n_elements(gstr)

  IF tag_exist(gstr[0], 'mu') THEN BEGIN 
      norms = gstr.p
      means = gstr.mu
      covs = gstr.cov
  ENDIF ELSE BEGIN 
      norms = gstr.norm
      means = gstr.mean
      covs = gstr.cov
  ENDELSE 

  xi = arrscl( findgen(nx), minx, maxx )

  model = fltarr(nx)
  
  FOR i=0L, ngauss-1 DO BEGIN 
      tmodel = norms[i]*gaussprob(xi, means[i], sqrt(covs[i]))
      
      model[*] = model[*] + tmodel[*]
  ENDFOR 

  norm = qgauss(model, xi, 200)

  
  model = model/norm
  model_struct = {ngauss:ngauss, $
                  x: xi, $
                  model: model}
  FOR i=0L, ngauss-1 DO BEGIN 

      name = 'model'+strn(i+1)

      tmodel = norms[i]*gaussprob(xi, means[i], sqrt(covs[i]))
      tmodel = tmodel/norm

      model_struct = $
        create_struct(model_struct, $
                      name, tmodel)

  ENDFOR 

  return,model_struct

END 


FUNCTION em::glike_1d, x, mu, cov, chi2min
;evaluate the likelihood at point x
;of a gaussian with mean mu and covariance cov
;x may be one dimensional only

  if n_params() LT 3 then begin
      message,'-syntax like = em->glike_1d(x, mu, cov, chi2min)'
  endif
  
  invcov=1.0/cov
  
  s=size(x) 
  
  chi2=(x-mu)*invcov*(x-mu)
  chi2=chi2-chi2min
  chi2=chi2 < 40.0
  like=exp(-.5*chi2)/(sqrt(cov)*(2*!pi)^2)
  
  return,like
end



FUNCTION em::gmix_1d, vect, ngauss, niter, $
           guess=guess, weight=weight, doplot=doplot, bin=bin

  if n_params() eq 0 then begin
      message, '-syntax gstr = em->gmix(vect,ngauss,niter,guess=, weight=, /doplot, bin=)'
  endif

  s=size(vect)
  ndim=1
  npoints=s(1)

  if n_elements(weight) ne 0 then begin
      ww=weight
  endif else begin
      ww=replicate(1.0/npoints,npoints)
  endelse
  weight=ww

  IF n_elements(guess) NE 0 THEN BEGIN 
      gstr = guess
  ENDIF ELSE BEGIN 

      gstr={norm:1.0,mean:0.0,cov:0.0}
      gstr=replicate(gstr,ngauss)
      mx=max(vect)
      mn=min(vect)
      gstr.mean=randomu(seed,ngauss)*(mx-mn)+mn
      gstr.cov=replicate((mx-mn)/((5.0*ngauss)^2),ngauss)+randomu(seed,ngauss)
  
  ENDELSE 

  if ngauss eq 1 then begin
      stats=moment(vect)
      gstr.mean=stats(0)
      gstr.cov=stats(1)
      gstr.norm=total(vect)
      return, gstr
  endif

  
  tau=fltarr(npoints,ngauss)
  beta=tau

  chi2m=replicate(100.0,npoints)
  for iter=0, niter-1 do begin
      
      for j=0,ngauss-1 do begin       
          like = self->glike_1d(vect, gstr[j].mean, gstr[j].cov, chi2m)
          beta[*,j]=like
      endfor  

      toti=beta#gstr.norm
      amat=(1.0/toti)#gstr.norm
      tau=beta*amat

      for j=0,ngauss-1 do begin
          gstr[j].norm=total(tau[*,j]*ww)
          gstr[j].mean=total(tau[*,j]*ww*vect)/gstr[j].norm
          gstr[j].cov=total(tau[*,j]*ww*(vect-gstr[j].mean)^2)/gstr[j].norm
      endfor                  
      
      IF ( ((iter+1) MOD 2) EQ 0 AND $
           keyword_set(doplot) AND $
           iter GT 0 ) THEN BEGIN 
          plothist, vect, xhist, yhist, bin=bin, /norm, $
            title = 'Iter = '+ntostr(iter+1)
          model_struct = self->gmix_model_1d(gstr, mn, mx, 1000)
          self->gmix_plot_model_1d,model_struct, /overplot

      ENDIF 

  endfor

  return, gstr

end



























FUNCTION em::cleanup

  return,1
  
END 


PRO em__define

  struct = { em, $
             dummy: 0 $
           }

END 
