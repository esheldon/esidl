      subroutine ad_mom(argc, argv)

c     This wrapper program is necessary to call the function from IDL
c     Compile with ifort -shared -static-intel -o ad_momi.so ad_momi.f90
c     Note the symbol will probably be ad_mom_

      integer*8 argc
      integer*8 argv(*)

      real*4 ad_mom2

      call ad_mom2( %val(argv(1)),%val(argv(2)),%val(argv(3)),
     & %val(argv(4)),%val(argv(5)),%val(argv(6)),%val(argv(7)),
     & %val(argv(8)),%val(argv(9)),%val(argv(10)),%val(argv(11)),
     & %val(argv(12)),%val(argv(13)),%val(argv(14)),
     & %val(argv(15)),%val(argv(16)),%val(argv(17)),
     & %val(argv(18)), %val(argv(19)))

      return
      end



      subroutine ad_mom2(ax,ay,ixx,iyy,ixy,nel,shiftmax,data,
     &sky,nx,ny,uncer,sigsky,mag,numiter,wcenx,wceny,whyflag,rho4)

c     calculates adaptive moments
c     interpolates weights for poorly sampled objects

      real tol1,tol2,xinterp2,xl,xh,yl,yh,xxx,yyy,xy,td,e1,e2,xinterp
      integer maxit

      
      integer nel,nx,ny
      real ax(nel),ay(nel),mag(nel)
      real ixx(nel),iyy(nel),ixy(nel)
      real sky(nel),uncer(nel),sigsky(nel),shiftmax,rho4(nel)
      integer imom
      integer ix1,ix2,iy1,iy2,i,j
      real w(2,2),m(2,2),n(2,2) ! ,sigxold,sigyold
      real xcen,ycen
      real sumx,sumy,grad,expon,weight,detm,detw,detn
      real xcenorig,ycenorig,xx1,yy1,sums4
      real data(nx,ny),spi,sumxx,sumyy,sumxy,sum,w1,w2,w12,xx,xx2,yy,yy2
      integer numiter(nel),whyflag(nel),kk
      real wcenx(nel),wceny(nel),e1old,e2old,ymod,m11old
      logical interpflag
      parameter(tol1=0.01,tol2=0.01,maxit=100,xinterp=3.0,
     &xinterp2=xinterp*xinterp)

c     write some info so you know it is here
c     and variables have been passed correctly
      
c     The main do loop now

c      write(*,*)shiftmax

      spi=sqrt(atan(1.)*4.)
      do kk=1,nel
        w(1,1)=max(1.,ixx(kk))
        w(2,2)=max(1.,iyy(kk))
        w(1,2)=ixy(kk)
c     the starting moments of the weight
        e1old=10.
        e2old=10.
        m11old=1.e6
        xcen=ax(kk)
        ycen=ay(kk)
        xcenorig=xcen
        ycenorig=ycen
c     the starting centroids
        
        imom=0
        do while(imom.lt.maxit)
          interpflag=.false.
          imom=imom+1
c     the number of iterations
          grad=4.*sqrt(max(w(1,1),w(2,2)))
          ix1=nint(max(xcen-grad-0.5,1.))
          iy1=nint(max(ycen-grad-0.5,1.))
          ix2=nint(min(xcen+grad+0.5,float(nx)))
          iy2=nint(min(ycen+grad+0.5,float(ny)))
c     defined the neighborhood of the object
          sumxx=0.
          sumyy=0.
          sumxy=0.
          sumx=0
          sumy=0.
          sum=0.
          sums4=0.
c          rmax=0.
c     loop over the neighborhood pixels
          detw=w(1,1)*w(2,2)-w(1,2)*w(1,2)
          if(detw.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=9
            goto 9090
          endif
          if(w(1,1).lt.xinterp.or.w(2,2).lt.xinterp.or.
     &    detw.lt.xinterp2)interpflag=.true.

c          if(interpflag)then
c            write(*,*)'here'
c          else
c            write(*,*)'there'
c          endif
          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.375
            xh=xx+0.375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.375
                yh=yy+0.375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then
                  do xxx=xl,xh,0.2499
                    xx2=xxx*xxx
                    do yyy=yl,yh,0.2499
                      yy2=yyy*yyy
                      expon=(xx2*w2+yy2*w1-2.*xxx*yyy*w12)
                      weight=exp(-0.5*expon)
                      ymod=(data(i,j)-sky(kk))*weight/16.
                      sumx=sumx+ymod*float(i)
                      sumy=sumy+ymod*float(j)
                      sum=sum+ymod
                    enddo
                  enddo
                endif
              else
                expon=xx2*w2+yy2*w1-2.*xx*yy*w12
                if(expon.le.9..and.expon.gt.0.)then
c     rmax=max(rmax,xx2+yy2)
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumx=sumx+ymod*float(i)
                  sumy=sumy+ymod*float(j)
                  sum=sum+ymod
                endif
              endif
            enddo
          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=1
            goto 9090
          endif
          xcen=sumx/sum
          ycen=sumy/sum
          if(abs(xcen-xcenorig).gt.shiftmax.or.
     &    abs(ycen-ycenorig).gt.shiftmax)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2
            goto 9090
          endif
          
          sum=0.
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.375
            xh=xx+0.375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.375
                yh=yy+0.375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then
                  do xxx=xl,xh,0.2499
                    xx2=xxx*xxx
                    do yyy=yl,yh,0.2499
                      yy2=yyy*yyy
                      expon=(xx2*w2+yy2*w1-2.*xxx*yyy*w12)
                      weight=exp(-0.5*expon)/16.
                      ymod=(data(i,j)-sky(kk))*weight
                      sumxx=sumxx+xx2*ymod
                      sumyy=sumyy+yy2*ymod
                      sumxy=sumxy+xxx*yyy*ymod
                      sums4=sums4+expon*expon*ymod
                      sum=sum+ymod
                    enddo
                  enddo              
                endif
              else
                xy=xx*yy
                expon=xx2*w2+yy2*w1-2.*xy*w12
                if(expon.le.9..and.expon.gt.0.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumxx=sumxx+xx2*ymod
                  sumyy=sumyy+yy2*ymod
                  sumxy=sumxy+xy*ymod
                  sums4=sums4+expon*expon*ymod
                  sum=sum+ymod
                endif
              endif
            enddo
          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=3
            goto 9090
          endif
          
          m(1,1)=sumxx/sum
          m(2,2)=sumyy/sum
          m(1,2)=sumxy/sum
          if(m(1,1).le.0..and.m(2,2).le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=4
            goto 9090
          endif
          td=w(1,1)+w(2,2)
          e1=(w(1,1)-w(2,2))/td
          e2=2.*w(1,2)/td
          if(abs(e1-e1old).lt.tol1.and.abs(e2-e2old).lt.tol1.and.
     &    abs(m(1,1)/m11old-1.).le.tol2)then
            ixx(kk)=w(1,1)
            iyy(kk)=w(2,2)
            ixy(kk)=w(1,2)
            rho4(kk)=sums4/sum
            ax(kk)=xcen
            ay(kk)=ycen
            detw=((w(1,1)*w(2,2)-w(1,2)*w(1,2)))**0.25
            whyflag(1)=0
            if(4.*sum-sums4.gt.0.)then
              uncer(kk)=4.*spi*sigsky(kk)*detw/(4.*sum-sums4)
            else
              uncer(kk)=9.999
            endif
            wcenx(kk)=xcen
            wceny(kk)=ycen
            numiter(kk)=imom
            goto 9090
          else
            detm=(m(1,1)*m(2,2)-m(1,2)*m(1,2))
c     write(*,*)'detm = ',detm,e1,e2,sumxx,sumyy,sumxy
            if(detm.le.1.e-7)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=5
              goto 9090
            endif
            
            detm=1./detm
            detw=1./detw
            n(1,1)=m(2,2)*detm-w(2,2)*detw
            n(2,2)=m(1,1)*detm-w(1,1)*detw
            n(1,2)=-m(1,2)*detm+w(1,2)*detw
            detn=n(1,1)*n(2,2)-n(1,2)*n(1,2)
            
c     write(*,*)'detn = ',detn,detw,
            if(detn.le.0.)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=6
              goto 9090
            endif
            
            detn=1./detn
            w(1,1)=n(2,2)*detn
            w(2,2)=n(1,1)*detn
            w(1,2)=-n(1,2)*detn
            e1old=e1
            e2old=e2
            m11old=m(1,1)
          endif
          if(w(1,1).lt.0..or.w(2,2).lt.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=7
            goto 9090
          endif
        enddo
        
        if(imom.eq.maxit)then
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
          numiter(kk)=imom
          whyflag(kk)=8
          goto 9090
        endif
        
 9090 enddo
      
      return
      
      end
      
      subroutine setbad(ixx,iyy,ixy,rho4,uncer)
      
      real ixx,iyy,ixy,rho4,uncer

      ixx=0.
      iyy=0.
      ixy=0.
      rho4=0.
      uncer=9.999

      return

      end
