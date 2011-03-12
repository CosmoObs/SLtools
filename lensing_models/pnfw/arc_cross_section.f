      program arc_cross_section
	
      double precision ks,e, lw_u, v_in,muth
c       ks: characteristic convergence of the PNFW 
c       e: parameter varepsilon of the PNFW
c       lw_u: length-to-width ratio threshold
c	ivt= input value for the  calculations
c	muth= threshold magnification (to use in the calculation of the sigma with the constraint u>u_th)
      double precision qu_min,qu_max,d_qu,qu
c     qu_min: minimum value of the quantity to be varied
c     qu_max: maximum value fo the quantity to be varied
c     d_qu: the step for the variation
      double precision v_out, sigma_q, sigma_0,sigma_u,u_min,sigma_u0
c       v_out: output value, only for gain speed
c	sigma_q: dimensionless arc cross section  (is necessary use iflag=2)
c	sigma_u: dimensionless arc cross section with the constraint u>u_th (is necessary use iflag=3)
c	u_min: minimum magnification
      double precision tinit, tend
c	tinit and tend only are useful for take control of the time.
      integer iflag,npt,iflag2,i,i_qu,n_qu,eflag
c	if iflag=1, this program calculate the arc cross section varying mu_th
c	if iflag=2, this program calculate the arc cross section varying ks
c	if iflag=3, this program calculate the arc cross section varying varepsilon
c	if iflag=4, this program calculate the arc cross section varying the length-to-width ratio	
c	iflag is useful to compute the arc cross section with/withour constrain on mu_th
c c
      double precision x1t,x2,ivt
c	x1t: Intersection of the tangential critical curve with the x1 axis at varepsilon =0
!       double precision sigma0,sigmae,kc,kcr,kcrit
!       integer ifil, ic
!       dimension sigma0(10000),sigmae(10000),kc(10000)

      open(unit=1,file='ifile',status='unknown')
      open(unit=2,file='ifile2',status='unknown')
      open(unit=100,file='tempo_calculo.txt',status='unknown',
     &    access='append')

      read(1,*)iflag, ks, e, lw_u, muth,eflag
!       read(1,*)iflag, ks, i, lw_u, muth !uncomment for the general case
!        e=(i-1)*0.025d0 ! uncomment in the general case
      read(2,*)qu_min, qu_max,n_qu

      d_qu=(qu_max-qu_min)/dfloat(n_qu)
c


      write(*,*)iflag, ks, e, lw_u, muth
      write(*,*)qu_min, qu_max, d_qu
!       pause
          
      write(*,*)"defining the value of the input value"

      if(iflag.eq.2)ks=qu_min
      if(ks.le.0.075d0)ivt=0.005d0
      if(ks.gt.0.075.and.ks.le.0.11d0)ivt=0.05d0
      if(ks.gt.0.11d0.and.ks.le.0.5d0)ivt=0.1d0
      if(ks.gt.0.5d0.and.ks.le.1.d0)ivt=0.50
      if(ks.gt.1.d0)ivt=1.0d0	
      write(*,*)ivt
c
c
        tinit=0.0d0
        call cpu_time(tinit)
        npt=251
c        
        
        x2=0.d0
        call raiz_lt(1,ivt,ks,0.d0,x1t,x2,eflag)
        ivt=x1t
        v_in=ivt
        iflag2=2
        do i_qu=1,n_qu

	  qu=qu_min +(i_qu-1)*d_qu
	  print*,i_qu,qu

          if(iflag.eq.1)then
            muth=qu
            iflag2=3
          endif
c
          if(iflag.eq.2)then
              ks=qu
              if(ks.le.0.075d0)ivt=0.005d0
              if(ks.gt.0.075.and.ks.le.0.11d0)ivt=0.05d0
              if(ks.gt.0.11d0.and.ks.le.0.5d0)ivt=0.1d0
              if(ks.gt.0.5d0.and.ks.le.1.d0)ivt=0.50
              if(ks.gt.1.d0)ivt=1.5d0	
           endif 
c
          if(iflag.eq.3)then
              call c_section(npt,v_in,ks,0.d0,lw_u,muth,iflag2,
     &	v_out,sigma_0,sigma_u0,u_min,eflag)
            v_in=v_out
              e=qu
          endif

          if(iflag.eq.4)lw_u=qu
c 
          
          call c_section(npt,v_in,ks,e,lw_u,muth,iflag2,
     &	v_out,sigma_q,sigma_u,u_min,eflag)
            v_in=v_out
            write(*,*)e,ks,lw_u,muth,sigma_q,sigma_u,u_min
            if(iflag.eq.1)write(3,*)muth,sigma_u/sigma_q
            if(iflag.eq.2)write(3,*)ks,sigma_q,u_min
!             if(iflag.eq.3)write(3,*)e,sigma_q/sigma_0,u_min
            if(iflag.eq.3)write(3,*)e,sigma_q,u_min
            if(iflag.eq.4)write(3,*)lw_u,sigma_q,u_min
!             write(4,*)e,ks,sigma_q
!             write(*,*)ks,e,sigma_q
        end do
        write(*,*)'The output files are fort.3'


	
!       if(iflag.eq.4)then
!         write(*,*)' calculating the arc cross section as function of ks'
!         tinit=0.0d0
!         call cpu_time(tinit)
!         npt=251
!         write(*,*)'Input the length-to-width ratio threshold'
!         read(*,*)lw_u
!         write(*,*)'Input the parameter varepsilon'
!         read(*,*) e	
!         dks=1.d-2
!         do ks=0.05d0,2.01d0,dks
!           write(*,*)'input value',ivt
!           el=e_arcs_lower(ks)-0.002
!           eu=e_arcs_upper(ks)+0.002
!           if(e.gt.el.and.e.lt.eu)then
!             sigma_q=-1.d0
!             write(*,*)'the arc cross section is not calculated'	
!           else
!             x2=0.d0
!             call raiz_lt(1,ivt,ks,0.d0,x1t,x2)
!             v_in=x1t
!             call c_section(npt,v_in,ks,e,lw_u,muth,2,
!     &	v_out,sigma_q,sigma_u,u_min)
!             ivt=v_out
!             write(*,*)ks,sigma_q,u_min
!             write(3,*)ks,sigma_q
!           endif
!         end do
!         write(*,*)'The output file is fort.3'
!       endif
! 
! 
!       if(iflag.eq.5)then
!         write(*,*)' calculating the arc cross section as function of varepsilon'
!         tinit=0.0d0
!         call cpu_time(tinit)
!         npt=251
!         write(*,*)'Input the length-to-width ratio threshold'
!         read(*,*)lw_u
!         write(*,*)'Input the characteristic convergence'
!         read(*,*) ks	
!         el=e_arcs_lower(ks)-0.002
!         eu=e_arcs_upper(ks)+0.02
!         if(ks.le.0.075d0)ivt=0.005d0
!         if(ks.gt.0.075.and.ks.le.0.11d0)ivt=0.05d0
!         if(ks.gt.0.11d0.and.ks.le.0.5d0)ivt=0.1d0
!         if(ks.gt.0.5d0.and.ks.le.1.d0)ivt=0.50
!         if(ks.gt.1.d0)ivt=1.5d0	
!         write(*,*)'the input value is',ivt
! ! 	  pause
!         de=0.001
!         x2=0.d0
!         call raiz_lt(1,ivt,ks,0.d0,x1t,x2)
!         ivt=x1t
!         v_in=ivt
!         do e=0.0d0,0.999d0,de
!           write(*,*)'input value',ivt
!           if(e.gt.el.and.e.lt.eu)then
!             sigma_q=-1.d0
!             write(*,*)'the arc cross section is not calculated'	
!           else
!             call c_section(npt,v_in,ks,e,lw_u,muth,2,
!     &	v_out,sigma_q,sigma_u,u_min)
!             v_in=v_out
!             write(*,*)ks,e,sigma_q,u_min
!           endif
!           write(4,*)ks,e,sigma_q
!         end do
!         write(*,*)'The output file is fort.4'
!       endif
! 
! 
!       if(iflag.eq.6)then
!         write(*,*)' calculating the arc cross section as function of lenght-to-width ratio'
!         tinit=0.0d0
!         call cpu_time(tinit)
!         npt=251
!         write(*,*)'Input the characteristic convergence'
!         read(*,*) ks	
!         el=e_arcs_lower(ks)-0.002
!         eu=e_arcs_upper(ks)+0.02
!         print*,el
!         write(*,*)'Input the parameter varepsilon'
!         read(*,*)e
!         dlw_u=0.1
!         write(*,*)'the input value is',ivt
!         pause
!         x2=0.d0
!         call raiz_lt(1,ivt,ks,0.d0,x1t,x2)
!         ivt=x1t
!         v_in=ivt
!         do lw_u=5.0d0,25.d0,dlw_u
!           write(*,*)'input value',ivt
!           if(e.gt.el.and.e.lt.eu)then
!             sigma_q=-1.d0
!             write(*,*)'the arc cross section is not calculated'	
!           else
!             call c_section(npt,v_in,ks,e,lw_u,muth,2,
!     &	v_out,sigma_q,sigma_u,u_min)
!             v_in=v_out
!             write(*,*)ks,e,sigma_q,u_min
!             write(5,*)lw_u,sigma_q
!           endif
!         end do
!         write(*,*)'The output file is fort.5'
!       endif

10    call cpu_time(tend)
      write(*,*)' Tempo :',tend-tinit,',seg'
      write(100,*)
      write(100,*)'ellipticity :',e
      write(100,*),'tempo pnfw :',tend-tinit,', seg'  
!       pause

      end