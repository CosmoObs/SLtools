! This file containt the intersted curves in the 
! lens plane, i.e critical curves and length-to-width
! ratio curves. And also make this plot in the source
! plane
        program contours
!         subroutine contours_pnfw(npt,lw,e,ks)
        double precision e,ks,lw,lwn,lwp,rs,theta
        double precision v_in,v_out,vit,vir,vilwn,vilwp
        double precision elp!,e_ad
        integer npt,opt,iflag
!         character*20 name
!         dimension name(8)
        open(unit=1,file='in_pnfw_par.txt',status='unknown')
!         open(unit=2,file='cluster_data.txt',status='unknown'
!      &  ,access='append')

        opt=2
        iflag=1
        lw=1.0
        read(1,*)ks,rs,elp,npt
!         read(1,*) name(1),ks
!         read(1,*) name(2),opt
!         read(1,*) name(3),elp
!         read(1,*) name(4),lw
!         read(1,*) name(5),npt
!         read(1,*) name(6),iflag
!         read(1,*) name(7),rs
        
!         
        if(opt.eq.2)then
            write(2,*)"working with the ellipticity parameter"
            e=elp
            write(*,*)e
        else
!             write(2,*)"working with the lens ellipticity"
!             e=e_ad(elp)
!             write(2,*)'the ellipticity of the angle deflection is', e
        endif
!         write(2,11)e
! 11      Format(5x,"ellipticity:",2x,1f6.3)
        write(*,*)"characteristic convergence",ks
!         write(*,*)name(2),opt
        write(*,*)"ellipticity",e
!         write(*,*)name(4),lw
        write(*,*)"scale radius",rs
        write(*,*)"number of points",npt
!         write(*,*)name(6),iflag
        
!         pause
!       Some scaling parameter
!         rs=1.d0
        theta=0.d0
!
        lwp=dabs(lw)
        lwn=-dabs(lw)
!
        if(ks.le.0.11d0)v_in=5.d-3
        if(ks.gt.0.11d0.and.ks.le.0.5d0)v_in=2.5d-1
        if(ks.gt.0.5d0.and.ks.le.1.0)v_in=0.5d0
        if(ks.gt.1.0)v_in=1.0d0

        vit=v_in
        
        call contour_lt(npt,vit,v_out,rs,theta,ks,e)
        if(iflag.eq.2)then
          vilwp=v_out
          vilwn=v_out
!           call contour_lwp(npt,vilwp,lwp,rs,theta,ks,e)  

!           call contour_lwn(npt,vilwn,lwn,rs,theta,ks,e,v_out)  
          go to 21
        endif
          vir=v_out
!         call contour_lr(npt,vir,rs,theta,ks,e)
!        
21      write(*,*)'Ending the plots'	
        end

****************************************************************************
****************************************************************************
        subroutine contour_lt(npt,v_in,v_out,rs,theta,ks,e)
        double precision v_in,v_out,ks,e
        double precision e0,x1,x2
        double precision vinter
        double precision x1t,x2t,x1tmax,stept,vit
        double precision y1t,y2t,y1,y2
        integer npt,i,eixo
        dimension x1t(1000),x2t(1000),y1t(1000),y2t(1000)
!       Useful to obtain the furthermore points
        double precision x2max1,x2max,x1ft,x2ft,stepft,v_ft
        integer i1,j,ntot
!       Useful for scaling coordinates
        double precision x1sc,x2sc,y1sc,y2sc,rs,theta,pi,thetar
        double precision rx1,rx2
        dimension x1sc(1000),x2sc(1000),y1sc(1000),y2sc(1000)
!       Useful for the physical limit
        double precision x1p,x2p,d2s,d2i,x1b,gamma2_e,y1p,y2p
        integer ip
        dimension x1p(1000),x2p(1000),y1p(1000),y2p(1000)
        external gamma2_e
      
!        
        pi=3.14159265d0
!       To convert from sexagesimal to radian angles		        
        thetar=theta*pi/180.d0
        rx1=rs*dcos(thetar)
        rx2=rs*dsin(thetar)
        
        eixo=1
        x2=0.d0
        e0=0.d0
        call raiz_lt(eixo,v_in,ks,e0,x1,x2)
        v_in=x1
        v_out=x1
!        
        eixo=1
        x2=0.d0
        vit=v_in
        call raiz_lt(eixo,vit,ks,e,x1,x2)
        x1tmax=x1
        x1p(npt)=x1tmax
        x2p(npt)=0.d0
        vinter=x1tmax
        call s_plane(ks,e,x1p(npt),x2p(npt),y1,y2)
        y1p(npt)=y1
        y2p(npt)=y2
        stept=x1p(npt)/dfloat(npt-1)
!	
        eixo=2
        vit=0.d0
        do i=npt-1,1,-1
          j=i-1
          x1=0.d0+(i-1)*stept
          call raiz_lt(eixo,vit,ks,e,x1,x2)
          x1p(i)=x1
          x2p(i)=x2
          vit=x2
          call s_plane(ks,e,x1p(i),x2p(i),y1,y2)
          y1p(i)=y1
          y2p(i)=y2
!
          x1b=0.d0+(j-1)*stept
          d2s=1.1*x2
          d2i=0.9*x2
!           if(gamma2_e(e,x1b,d2s).lt.0.d0.or.
!      &    gamma2_e(e,x1b,d2i).lt.0.d0)then
!             write(*,*)"limit parameters in lt= 0", ks,e
!             write(60,*)'tangential critical curve ',ks,e
!             go to 10
!           endif
        end do 
10      ip=i
        if(ip.lt.1)ip=1
!         write(*,*)'printing the lower value',ip
!      Redefining the number of points
        do i=npt,ip,-1
            j=i-ip+1
          x1t(j)=x1p(i)
          x2t(j)=x2p(i)
          y1t(j)=y1p(i)
          y2t(j)=y2p(i)
        end do
!
        npt=npt-ip+1


        x2max1=x2t(1)
        do i=1,npt
            i1=i+1
            if(x2t(i1).gt.x2t(i))then
              x2max1=x2t(i1)
            endif
        end do
        x2max=x2max1
!
        j=npt
        stepft=x2max/dfloat(npt)
        i1=1
        v_ft=vinter
20      x2ft=i1*stepft
        eixo=1
        call raiz_lt(eixo,v_ft,ks,e,x1ft,x2ft)
        if(x1ft.ge.0.99d0*x1t(npt))then
            v_ft=x1ft
            if(x1ft.gt.x1t(npt))then
	j=j+1
	x1t(j)=x1ft
	x2t(j)=x2ft
	call s_plane(ks,e,x1t(j),x2t(j),y1,y2)
	y1t(j)=y1
	y2t(j)=y2
            endif
            i1=i1+1  
            go to 20
        endif
        ntot=j

!       Scaling our lengths
        do i=ntot,1,-1
          x1sc(i)=x1t(i)*rx1+x2t(i)*rx2
          x2sc(i)=-x1t(i)*rx2+x2t(i)*rx1
          y1sc(i)=y1t(i)*rx1+y2t(i)*rx2
          y2sc(i)=-y1t(i)*rx2+y2t(i)*rx1
        end do

!       Writing in the first quadrant
        write(61,*)x1sc(npt),x2sc(npt)
        write(65,*)x1sc(npt),x2sc(npt)
        write(71,*)y1sc(npt),y2sc(npt)
        write(75,*)y1sc(npt),y2sc(npt)
        do i=npt+1,ntot
            write(61,*)x1sc(i),x2sc(i)
            write(65,*)x1sc(i),x2sc(i)
            write(71,*)y1sc(i),y2sc(i)
            write(75,*)y1sc(i),y2sc(i)
        end do

        do i=npt-1,1,-1
            write(61,*)x1sc(i),x2sc(i)
            write(65,*)x1sc(i),x2sc(i)
            write(71,*)y1sc(i),y2sc(i)
            write(75,*)y1sc(i),y2sc(i)
        end do
!
!      Writing in the second quadrant        
        do i=1,npt-1
            write(62,*)-x1sc(i),x2sc(i)
            write(65,*)-x1sc(i),x2sc(i)
            write(72,*)-y1sc(i),y2sc(i)
            write(75,*)-y1sc(i),y2sc(i)

        end do

        do i=ntot,npt,-1
            write(62,*)-x1sc(i),x2sc(i)
            write(65,*)-x1sc(i),x2sc(i)
            write(72,*)-y1sc(i),y2sc(i)
            write(75,*)-y1sc(i),y2sc(i)
        end do

!       Writing in the thirt quadrant
        write(63,*)-x1sc(npt),-x2sc(npt)
        write(65,*)-x1sc(npt),-x2sc(npt)
        write(73,*)-y1sc(npt),-y2sc(npt)
        write(75,*)-y1sc(npt),-y2sc(npt)

        do i=npt+1,ntot
            write(63,*)-x1sc(i),-x2sc(i)
            write(65,*)-x1sc(i),-x2sc(i)
            write(73,*)-y1sc(i),-y2sc(i)
            write(75,*)-y1sc(i),-y2sc(i)
        end do

        do i=npt-1,1,-1
            write(63,*)-x1sc(i),-x2sc(i)
            write(65,*)-x1sc(i),-x2sc(i)
            write(73,*)-y1sc(i),-y2sc(i)
            write(75,*)-y1sc(i),-y2sc(i)
        end do
!        Writing in the fourt quadrant
        do i=1,npt-1
            write(64,*)x1sc(i),-x2sc(i)
            write(65,*)x1sc(i),-x2sc(i)
            write(74,*)y1sc(i),-y2sc(i)
            write(75,*)y1sc(i),-y2sc(i)
        end do

        do i=ntot,npt,-1
            write(64,*)x1sc(i),-x2sc(i)
            write(65,*)x1sc(i),-x2sc(i)
            write(74,*)y1sc(i),-y2sc(i)
            write(75,*)y1sc(i),-y2sc(i)
        end do  

        return
        
        end
! ****************************************************************************
! ****************************************************************************
        subroutine contour_lwn(npt,v_in,lwn,rs,theta,ks,e,v_out)
        double precision v_in,lwn,ks,e,v_out
        double precision e0,x1,x2
        double precision x1lw,x2lw,x1lwmax,steplw,vilw
        double precision y1lw,y2lw,y1,y2
        integer i,npt,eixo
        dimension x1lw(1000),x2lw(1000),y1lw(1000),y2lw(1000)
!       Useful to obtain the furthermore points
        double precision x2max1,x2max,x1ft,x2ft,stepft,v_ft,vinter
        integer i1,j,ntot
!       Useful for scaling coordinates
        double precision x1sc,x2sc,y1sc,y2sc,rs,theta,pi,thetar
        double precision rx1,rx2
        dimension x1sc(1000),x2sc(1000),y1sc(1000),y2sc(1000)
!       Useful for the physical limit
        double precision x1p,x2p,d2s,d2i,x1b,gamma2_e,y1p,y2p
        integer ip
        dimension x1p(1000),x2p(1000),y1p(1000),y2p(1000)
        external gamma2_e
!        
        pi=3.14159265d0
!       To convert from sexagesimal to radian angles		        
        thetar=theta*pi/180.d0
        rx1=rs*dcos(thetar)
        rx2=rs*dsin(thetar)
!         write(*,*)'scaling coordinates',rx1,rx2


        eixo=1
        x2=0.d0
        e0=0.d0
        call raiz_lw(eixo,lwn,v_in,ks,e0,x1,x2) ! v_in is from v_out of the tangential critical curve
        v_in=x1
        v_out=x1 ! it's useful for calculate the radial critical curve
!
        vilw=v_in
        call raiz_lw(eixo,lwn,vilw,ks,e,x1,x2)
        x1lwmax=x1
        x1p(npt)=x1lwmax
        x2p(npt)=0.d0
        vinter=x1lwmax
!
        call s_plane(ks,e,x1p(npt),x2p(npt),y1,y2)
        y1p(npt)=y1
        y2p(npt)=y2
!
        steplw=x1p(npt)/dfloat(npt-1)
        vilw=0.d0
        eixo=2  
        do i=npt-1,1,-1
            x1=0.d0+(i-1)*steplw
            call raiz_lw(eixo,lwn,vilw,ks,e,x1,x2)
            x1p(i)=x1
            x2p(i)=x2
            vilw=x2
            call s_plane(ks,e,x1p(i),x2p(i),y1,y2)
            y1p(i)=y1
            y2p(i)=y2
!
            j=i-1
            x1b=0.d0+(j-1)*steplw
            d2s=1.05*x2
            d2i=0.95*x2
            if(gamma2_e(e,x1b,d2s).lt.0.d0.or.
     &      gamma2_e(e,x1b,d2i).lt.0.d0)then
              write(*,*)"parameters limit in L/W = -10", ks,e
              write(60,*)'L/W = -10  critical curve ',ks,e
              go to 10
            endif
        end do
10      ip=i
        if(ip.lt.1)ip=1
!         write(*,*)'printing the lower value',ip
!      Redefining the number of points
        do i=npt,ip,-1
          j=i-ip+1
          x1lw(j)=x1p(i)
          x2lw(j)=x2p(i)
          y1lw(j)=y1p(i)
          y2lw(j)=y2p(i)
        end do
        npt=npt-ip+1
!
        x2max1=x2lw(1)
        do i=1,npt
!             if(i.eq.1)print*, 'o valor inicial for lw = -10',x2max1
            i1=i+1
            if(x2lw(i1).gt.x2lw(i))then
              x2max1=x2lw(i1)
            endif
        end do
        x2max=x2max1
!         write(*,*)' the maximum value for lw = -10 is',x2max
        j=npt
        stepft=x2max/dfloat(npt)
        i1=1
        v_ft=vinter
20      x2ft=i1*stepft
        eixo=1
        call raiz_lw(eixo,lwn,v_ft,ks,e,x1ft,x2ft)
        if(x1ft.ge.0.99d0*x1lw(npt))then
            v_ft=x1ft
            if(x1ft.gt.x1lw(npt))then
	j=j+1
	x1lw(j)=x1ft
	x2lw(j)=x2ft
	call s_plane(ks,e,x1lw(j),x2lw(j),y1,y2)
	y1lw(j)=y1
	y2lw(j)=y2
            endif
            i1=i1+1  
            go to 20
        endif
        ntot=j
        print*,'the total number of lw = -10 is',ntot,npt
!         pause
!         write(*,*)
!         write(*,*)
!       Scaling our lengths
        do i=ntot,1,-1
          x1sc(i)=x1lw(i)*rx1+x2lw(i)*rx2
          x2sc(i)=-x1lw(i)*rx2+x2lw(i)*rx1
          y1sc(i)=y1lw(i)*rx1+y2lw(i)*rx2
          y2sc(i)=-y1lw(i)*rx2+y2lw(i)*rx1
        end do

!       Writing in the first quadrant
        write(41,*)x1sc(npt),x2sc(npt)
        write(45,*)x1sc(npt),x2sc(npt)
        write(51,*)y1sc(npt),y2sc(npt)
        write(55,*)y1sc(npt),y2sc(npt)
!         do i=npt,ntot
!             write(41,*)x1sc(i),x2sc(i)
!             write(45,*)x1sc(i),x2sc(i)
!             write(51,*)y1sc(i),y2sc(i)
!             write(55,*)y1sc(i),y2sc(i)
!         end do

        do i=npt-1,1,-1
            write(41,*)x1sc(i),x2sc(i)
            write(45,*)x1sc(i),x2sc(i)
            write(51,*)y1sc(i),y2sc(i)
            write(55,*)y1sc(i),y2sc(i)
        end do
!      Writing in the second quadrant        
        do i=1,npt-1
            write(42,*)-x1sc(i),x2sc(i)
            write(45,*)-x1sc(i),x2sc(i)
            write(52,*)-y1sc(i),y2sc(i)
            write(55,*)-y1sc(i),y2sc(i)

        end do

!         do i=ntot,npt,-1
!             write(42,*)-x1sc(i),x2sc(i)
!             write(45,*)-x1sc(i),x2sc(i)
!             write(52,*)-y1sc(i),y2sc(i)
!             write(55,*)-y1sc(i),y2sc(i)
!         end do

!       Writing in the thirt quadrant
        write(43,*)-x1sc(npt),-x2sc(npt)
        write(45,*)-x1sc(npt),-x2sc(npt)
        write(53,*)-y1sc(npt),-y2sc(npt)
        write(55,*)-y1sc(npt),-y2sc(npt)

!         do i=npt+1,ntot
!             write(43,*)-x1sc(i),-x2sc(i)
!             write(45,*)-x1sc(i),-x2sc(i)
!             write(53,*)-y1sc(i),-y2sc(i)
!             write(55,*)-y1sc(i),-y2sc(i)
!         end do

        do i=npt-1,1,-1
            write(43,*)-x1sc(i),-x2sc(i)
            write(45,*)-x1sc(i),-x2sc(i)
            write(53,*)-y1sc(i),-y2sc(i)
            write(55,*)-y1sc(i),-y2sc(i)
        end do
!        Writing in the fourt quadrant
        do i=1,npt
            write(44,*)x1sc(i),-x2sc(i)
            write(45,*)x1sc(i),-x2sc(i)
            write(54,*)y1sc(i),-y2sc(i)
            write(55,*)y1sc(i),-y2sc(i)
        end do

!         do i=ntot,npt,-1
!             write(44,*)x1sc(i),-x2sc(i)
!             write(45,*)x1sc(i),-x2sc(i)
!             write(54,*)y1sc(i),-y2sc(i)
!             write(45,*)y1sc(i),-y2sc(i)
!         end do  

        return
        
        end

****************************************************************************
****************************************************************************
        subroutine contour_lwp(npt,v_in,lwp,rs,theta,ks,e)
        double precision v_in,lwp,ks,e
        double precision e0,x1,x2
        double precision x1lw,x2lw,x1lwmax,steplw,vilw
        double precision y1lw,y2lw,y1,y2
        integer i,npt,eixo
        dimension x1lw(1000),x2lw(1000),y1lw(1000),y2lw(1000)
!       Useful to obtain the furthermore points
        double precision x2max1,x2max,x1ft,x2ft,stepft,v_ft,vinter
        integer i1,j,ntot
!       Useful for scaling coordinates
        double precision x1sc,x2sc,y1sc,y2sc,rs,theta,pi,thetar
        double precision rx1,rx2
        dimension x1sc(1000),x2sc(1000),y1sc(1000),y2sc(1000)
!       Useful for the physical limit
        double precision x1p,x2p,d2s,d2i,x1b,gamma2_e,y1p,y2p
        integer ip
        dimension x1p(1000),x2p(1000),y1p(1000),y2p(1000)
        external gamma2_e

!        
        pi=3.14159265d0
!       To convert from sexagesimal to radian angles		        
        thetar=theta*pi/180.d0
        rx1=rs*dcos(thetar)
        rx2=rs*dsin(thetar)
!         write(*,*)'scaling coordinates',rx1,rx2  


        npt=npt
        eixo=1
        x2=0.d0
        e0=0.d0
        call raiz_lw(eixo,lwp,v_in,ks,e0,x1,x2) ! v_in is from v_out of the tangential critical curve
        v_in=x1
!
        vilw=v_in
        call raiz_lw(eixo,lwp,vilw,ks,e,x1,x2)
        x1lwmax=x1
        x1p(npt)=x1lwmax
        x2p(npt)=0.d0
        vinter=x1lwmax
!
        call s_plane(ks,e,x1p(npt),x2p(npt),y1,y2)
        y1p(npt)=y1
        y2p(npt)=y2
!
        steplw=x1p(npt)/dfloat(npt-1)
        vilw=0.d0
        eixo=2  
        do i=npt-1,1,-1
            x1=0.d0+(i-1)*steplw
            call raiz_lw(eixo,lwp,vilw,ks,e,x1,x2)
            x1p(i)=x1
            x2p(i)=x2
            vilw=x2
            call s_plane(ks,e,x1p(i),x2p(i),y1,y2)
            y1p(i)=y1
            y2p(i)=y2
!
            j=i-1
            x1b=0.d0+(j-1)*steplw
            d2s=1.05*x2
            d2i=0.95*x2
            if(gamma2_e(e,x1b,d2s).lt.0.d0.or.
     &      gamma2_e(e,x1b,d2i).lt.0.d0)then
              write(*,*)"parameters limit in L/W = 10 ", ks,e
              write(60,*)'L/W = 10  critical curve ',ks,e
              go to 10
            endif
        end do
10      ip=i
!         write(*,*)'printing the lower value',ip
        if(ip.lt.1)ip=1
!      Redefining the number of points
        do i=npt,ip,-1
          j=i-ip+1
          x1lw(j)=x1p(i)
          x2lw(j)=x2p(i)
          y1lw(j)=y1p(i)
          y2lw(j)=y2p(i)
        end do
        npt=npt-ip+1
!
        x2max1=x2lw(1)
        do i=1,npt
!             if(i.eq.1)print*, 'o valor inicial for lw = 10',x2max1
            i1=i+1
            if(x2lw(i1).gt.x2lw(i))then
              x2max1=x2lw(i1)
            endif
        end do
        x2max=x2max1
!         write(*,*)' the maximum value for lw = 10 is',x2max
        j=npt
        stepft=x2max/dfloat(npt)
        i1=1
        v_ft=vinter
20      x2ft=i1*stepft
        eixo=1
        call raiz_lw(eixo,lwp,v_ft,ks,e,x1ft,x2ft)
        if(x1ft.ge.0.99d0*x1lw(npt))then
            v_ft=x1ft
            if(x1ft.gt.x1lw(npt))then
	j=j+1
	x1lw(j)=x1ft
	x2lw(j)=x2ft
	call s_plane(ks,e,x1lw(j),x2lw(j),y1,y2)
	y1lw(j)=y1
	y2lw(j)=y2
            endif
            i1=i1+1  
            go to 20
        endif
        ntot=j
!         print*,'the total number of points for lw = 10 is',ntot
!         write(*,*)
!         write(*,*)
!       Scaling our lengths
        do i=ntot,1,-1
          x1sc(i)=x1lw(i)*rx1+x2lw(i)*rx2
          x2sc(i)=-x1lw(i)*rx2+x2lw(i)*rx1
          y1sc(i)=y1lw(i)*rx1+y2lw(i)*rx2
          y2sc(i)=-y1lw(i)*rx2+y2lw(i)*rx1
        end do

!       Writing in the first quadrant
        write(81,*)x1sc(npt),x2sc(npt)
        write(85,*)x1sc(npt),x2sc(npt)
        write(91,*)y1sc(npt),y2sc(npt)
        write(95,*)y1sc(npt),y2sc(npt)
        do i=npt+1,ntot
            write(81,*)x1sc(i),x2sc(i)
            write(85,*)x1sc(i),x2sc(i)
            write(91,*)y1sc(i),y2sc(i)
            write(95,*)y1sc(i),y2sc(i)
        end do

        do i=npt-1,2,-1
            write(81,*)x1sc(i),x2sc(i)
            write(85,*)x1sc(i),x2sc(i)
            write(91,*)y1sc(i),y2sc(i)
            write(95,*)y1sc(i),y2sc(i)
        end do
!      Writing in the second quadrant        
        do i=1,npt-1
            write(82,*)-x1sc(i),x2sc(i)
            write(85,*)-x1sc(i),x2sc(i)
            write(92,*)-y1sc(i),y2sc(i)
            write(95,*)-y1sc(i),y2sc(i)

        end do

        do i=ntot,npt,-1
            write(82,*)-x1sc(i),x2sc(i)
            write(85,*)-x1sc(i),x2sc(i)
            write(92,*)-y1sc(i),y2sc(i)
            write(95,*)-y1sc(i),y2sc(i)
        end do

!       Writing in the thirt quadrant
        write(83,*)-x1sc(npt),-x2sc(npt)
        write(85,*)-x1sc(npt),-x2sc(npt)
        write(93,*)-y1sc(npt),-y2sc(npt)
        write(95,*)-y1sc(npt),-y2sc(npt)

        do i=npt+1,ntot
            write(83,*)-x1sc(i),-x2sc(i)
            write(85,*)-x1sc(i),-x2sc(i)
            write(93,*)-y1sc(i),-y2sc(i)
            write(95,*)-y1sc(i),-y2sc(i)
        end do

        do i=npt-1,1,-1
            write(83,*)-x1sc(i),-x2sc(i)
            write(85,*)-x1sc(i),-x2sc(i)
            write(93,*)-y1sc(i),-y2sc(i)
            write(95,*)-y1sc(i),-y2sc(i)
        end do
!        Writing in the fourt quadrant
        do i=1,npt-1
            write(84,*)x1sc(i),-x2sc(i)
            write(85,*)x1sc(i),-x2sc(i)
            write(94,*)y1sc(i),-y2sc(i)
            write(95,*)y1sc(i),-y2sc(i)
        end do

        do i=ntot,npt,-1
            write(84,*)x1sc(i),-x2sc(i)
            write(85,*)x1sc(i),-x2sc(i)
            write(94,*)y1sc(i),-y2sc(i)
            write(95,*)y1sc(i),-y2sc(i)
        end do  
!
        return
        
        end

****************************************************************************
****************************************************************************
        subroutine s_plane(ks,e,x1,x2,y1,y2)	
c      This subroutine make the mapping from the lens into source plane. 
c      It use the lens equation:\vec{y}=-vec{x}-\alpha(\vec{x})
        double precision ks,e,x1,x2,y1,y2
        double precision y1_s,y2_s
        external y1_s,y2_s
        y1=dabs(y1_s(ks,e,x1,x2))
        y2=dabs(y2_s(ks,e,x1,x2))
        return
        end	
****************************************************************************
****************************************************************************
        subroutine contour_lr(npt,v_in,rs,theta,ks,e)
        double precision v_in,ks,e
        double precision e0,x1,x2
        double precision x1r,x2r,x1rmax,stepr,vir
        double precision y1r,y2r,y1,y2
        integer npt,i,eixo
        dimension x1r(1000),x2r(1000),y1r(1000),y2r(1000)
!       Useful to obtain the furthermore points
        double precision x2max1,x2max,x1ft,x2ft,stepft,v_ft
        integer i1,j,ntot
!       Useful for scaling coordinates
        double precision x1sc,x2sc,y1sc,y2sc,rs,theta,pi,thetar
        double precision rx1,rx2
        dimension x1sc(1000),x2sc(1000),y1sc(1000),y2sc(1000)
!       Useful for the physical limit
        double precision x1p,x2p,d2s,d2i,x1b,gamma2_e,y1p,y2p
        integer ip
        dimension x1p(1000),x2p(1000),y1p(1000),y2p(1000)
!        
        pi=3.14159265d0
!       To convert from sexagesimal to radian angles		        
        thetar=theta*pi/180.d0
        rx1=rs*dcos(thetar)
        rx2=rs*dsin(thetar)
!
        eixo=1
        x2=0.d0
        e0=0.d0
        call raiz_lr(eixo,v_in,ks,e0,x1,x2)
!
        vir=x1
        call raiz_lr(eixo,vir,ks,e,x1,x2)
        x1rmax=x1
        x1p(npt)=x1rmax
        x2p(npt)=0.d0
        vinter=x1rmax  
        call s_plane(ks,e,x1p(npt),x2p(npt),y1,y2)
        y1p(npt)=y1
        y2p(npt)=y2
!        
        stepr=x1p(npt)/dfloat(npt-1)
        vir=0.d0
        eixo=2  
        do i=npt-1,1,-1
            j=i-1
            x1=0.d0+(i-1)*stepr
            call raiz_lr(eixo,vir,ks,e,x1,x2)
            x1p(i)=x1
            x2p(i)=x2
            vir=x2
            call s_plane(ks,e,x1p(i),x2p(i),y1,y2)
            y1p(i)=y1
            y2p(i)=y2
!
            x1b=0.d0+(j-1)*stepr
            d2s=1.1*x2
            d2i=0.9*x2
            if(gamma2_e(e,x1b,d2s).lt.0.d0.or.
     &      gamma2_e(e,x1b,d2i).lt.0.d0)then
              write(*,*)"parameters", ks,e
              write(60,*)'radial critical curve ',ks,e
              go to 10
            endif
        end do
10      ip=i
        if(ip.lt.1)ip=1
!         write(*,*)'printing the lower value',ip
!      Redefining the number of points
        do i=npt,ip,-1
          j=i-ip+1
          x1r(j)=x1p(i)
          x2r(j)=x2p(i)
          y1r(j)=y1p(i)
          y2r(j)=y2p(i)
        end do
        npt=npt-ip+1
!
        x2max1=x2r(1)
        do i=1,npt
            i1=i+1
            if(x2r(i1).gt.x2r(i))then
              x2max1=x2r(i1)
            endif
        end do
        x2max=x2max1
!
        j=npt
        stepft=x2max/dfloat(npt)
        i1=1
        v_ft=vinter
20      x2ft=i1*stepft
        eixo=1
        call raiz_lr(eixo,v_ft,ks,e,x1ft,x2ft)
        if(x1ft.ge.0.99d0*x1r(npt))then
            v_ft=x1ft
            if(x1ft.gt.x1r(npt))then
	j=j+1
	x1r(j)=x1ft
	x2r(j)=x2ft
	call s_plane(ks,e,x1r(j),x2r(j),y1,y2)
	y1r(j)=y1
	y2r(j)=y2
            endif
            i1=i1+1  
            go to 20
        endif
        ntot=j
!         print*,'the total number of points for lr = 0 is',ntot
!       Scaling our lengths
        do i=ntot,1,-1
          x1sc(i)=x1r(i)*rx1+x2r(i)*rx2
          x2sc(i)=-x1r(i)*rx2+x2r(i)*rx1
          y1sc(i)=y1r(i)*rx1+y2r(i)*rx2
          y2sc(i)=-y1r(i)*rx2+y2r(i)*rx1
        end do

!       Writing in the first quadrant
        write(21,*)x1sc(npt),x2sc(npt)
        write(25,*)x1sc(npt),x2sc(npt)
        write(31,*)y1sc(npt),y2sc(npt)
        write(35,*)y1sc(npt),y2sc(npt)
        do i=npt+1,ntot
            write(21,*)x1sc(i),x2sc(i)
            write(25,*)x1sc(i),x2sc(i)
            write(31,*)y1sc(i),y2sc(i)
            write(35,*)y1sc(i),y2sc(i)
        end do

        do i=npt-1,1,-1
            write(21,*)x1sc(i),x2sc(i)
            write(25,*)x1sc(i),x2sc(i)
            write(31,*)y1sc(i),y2sc(i)
            write(35,*)y1sc(i),y2sc(i)
        end do
!      Writing in the second quadrant        
        do i=1,npt-1
            write(22,*)-x1sc(i),x2sc(i)
            write(25,*)-x1sc(i),x2sc(i)
            write(32,*)-y1sc(i),y2sc(i)
            write(35,*)-y1sc(i),y2sc(i)

        end do

        do i=ntot,npt,-1
            write(22,*)-x1sc(i),x2sc(i)
            write(25,*)-x1sc(i),x2sc(i)
            write(32,*)-y1sc(i),y2sc(i)
            write(35,*)-y1sc(i),y2sc(i)
        end do

!       Writing in the thirt quadrant
        write(23,*)-x1sc(npt),-x2sc(npt)
        write(25,*)-x1sc(npt),-x2sc(npt)
        write(33,*)-y1sc(npt),-y2sc(npt)
        write(35,*)-y1sc(npt),-y2sc(npt)

        do i=npt+1,ntot
            write(23,*)-x1sc(i),-x2sc(i)
            write(25,*)-x1sc(i),-x2sc(i)
            write(33,*)-y1sc(i),-y2sc(i)
            write(35,*)-y1sc(i),-y2sc(i)
        end do

        do i=npt-1,1,-1
            write(23,*)-x1sc(i),-x2sc(i)
            write(25,*)-x1sc(i),-x2sc(i)
            write(33,*)-y1sc(i),-y2sc(i)
            write(35,*)-y1sc(i),-y2sc(i)
        end do
!        Writing in the fourt quadrant
        do i=1,npt-1
            write(24,*)x1sc(i),-x2sc(i)
            write(25,*)x1sc(i),-x2sc(i)
            write(34,*)y1sc(i),-y2sc(i)
            write(35,*)y1sc(i),-y2sc(i)
        end do

        do i=ntot,npt,-1
            write(24,*)x1sc(i),-x2sc(i)
            write(25,*)x1sc(i),-x2sc(i)
            write(34,*)y1sc(i),-y2sc(i)
            write(35,*)y1sc(i),-y2sc(i)
        end do  

        return
        
        end