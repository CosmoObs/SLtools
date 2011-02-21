!   THIS FILE CONTAINT ALL THE USEFUL FINDING-ROOT SUBROUTINES.
!   IN THEIR WE USE TWO METHOD*  HELPED WITH A SUBROUTINES THAT
!   FIND THE INTERVAL  WHERE THE POSSIBLE ROOT EXISTS.
!  (*) THE FIRST METHOD WE CALLED AS NEWTON-RAPSHON IN FINITE DIFFERENCES,
!   AND THE SECOND  IS THE WELL-KNOWN BISECTION METHOD.
!  ADVICE: THESE SUBROUTINE SHOULD BE VERIFIED BY OTHER SOFTWARE, i.e. MATHEMATICA

!   SUBROUTINE raiz_ke
!     This subroutine find the root of the convergence of the PNFW (following the angle 
!     deflection model), and use the following subroutines*
!	SUBROUTINE inter_ke : Find the range  
!	SUBROUTINE nrdf_ke : Find the root applying the firts method
!	SUBROUTINE bisec_ke: Find the root appliying the second ones
!   SUBROUTINE raiz_lr
!   This subroutine find the root of the radial eigenvalue of the matrix tensor,
!   and use the following subroutines*
!	SUBROUTINE inter_lr : Find the range  
!	SUBROUTINE nrdf_lr : Find the root appliying the firts method
!	SUBROUTINE bisec_lr: Find the root appliying the second ones


!   SUBROUTINE raiz_lt
!   This subroutine find the root of the tangential eigenvalue of the matrix tensor,
!   and use the following subroutines*
!	SUBROUTINE inter_lt : Find the range  
!	SUBROUTINE nrdf_lt : Find the root appliying the firts method
!	SUBROUTINE bisec_lt: Find the root appliying the second ones


!   SUBROUTINE raiz_lw
!   This subroutine find the root of the length-to-width ratio follows the ratio 
!   between the eigenvalues of the matrix tensor, and use the following subroutines*
!	SUBROUTINE inter_lw : Find the range  
!	SUBROUTINE nrdf_lw : Find the root appliying the firts method
!	SUBROUTINE bisec_lw: Find the root appliying the second ones

*************************************O O O*****************************************************

      SUBROUTINE  raiz_ke(eixo,v_input,iscont,e,x1,x2,eflag)
      double precision v_input,iscont,e,x1,x2
      double precision z_inf,z_sup
      integer eixo,iter,eflag
	
      call inter_ke(eixo,v_input,iscont,e,x1,x2,z_inf,z_sup,eflag)	
      call nrdf_ke(eixo,z_inf,z_sup,iscont,e,x1,x2,iter,eflag)
!      
      if(iter.gt.10)then
          write(*,*)'USING THE BISECTION METHOD'
          call bisec_ke(eixo,z_inf,z_sup,iscont,e,x1,x2,iter,eflag)
      endif
!
      x1=dabs(x1)
      x2=dabs(x2)	
      RETURN
      END
	
**********************************************************************************************
*********************************************************************************************	
      SUBROUTINE inter_ke(eixo,v_input,iscont,e,x1,x2,z_inf,z_sup,eflag)
!     This subroutine find the possible range in which the root can be exists
!     Input Data: v_input: Initial value	
!     Output Data: z_inf and z_sup are the  ()  and the upper limit of the range	
      double precision v_input,iscont,e,x1,x2,z_inf,z_sup
      double precision kappa_e2,z,dz,xi,xj
      integer eixo,eflag
      external kappa_e2
!	
      if(v_input.eq.0.d0)then
          z=1.d-3
      else
          z=v_input
      endif	
!
      dz=z*1.d-3
!	
      if(eixo.eq.1)then
          xi=z
          xj=x2
      else	
          xi=x1
          xj=z
      endif
!	
      if(kappa_e2(iscont,e,xi,xj,eflag).gt.0.d0)then
          z=z+dz
210     if(kappa_e2(iscont,e,xi,xj,eflag).gt.0.d0)then
            z_inf=z
            z=z+dz
            if(eixo.eq.1)then
	  xi=z
	  xj=x2
            else
	  xi=x1
	  xj=z
            endif	
!
          go to 210
!
        else
           z_sup=z
           if(eixo.eq.1)then
	xi=z_sup
	xj=x2
            else
	xi=x1
	xj=z
            endif
				
            go to 230

        endif	
      else	
          z=z-dz
220     if(kappa_e2(iscont,e,xi,xj,eflag).lt.0.d0)then
          z_sup=z
          z=z-dz
          if(eixo.eq.1)then
	xi=z
	xj=x2
          else
	xi=x1
	xj=z
          endif	
          go to 220
        else
          z_inf=z
          if(eixo.eq.1)then
	xi=z_inf
	xj=x2
          else
	xi=x1
	xj=z_inf
          endif				
          go to 230
        endif			
      endif

230   RETURN
      END
	
	
**********************************************************************
**********************************************************************
      SUBROUTINE nrdf_ke(eixo,z_inf,z_sup,iscont,e,x1,x2,i,eflag)
      double precision z_inf,z_sup,iscont,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj,zk,tol_in
      double precision kappa_e2,f0,f1,fm,gm,hm
      double precision del1,del2,del_min,fk
      INTEGER i,eixo,eflag
      EXTERNAL kappa_e2
c	
      i=0
!
      if(e.lt.1.d-1)then
          tol=1.d-4
          tol_in=1.d-4
      else
        tol=1.d-8
        tol_in=1.d-6
      endif	
!
      z0=z_inf
      z1=z_sup
!
      if(eixo.eq.1)then
        xj=x2
        xi=z0
        f0=kappa_e2(iscont,e,xi,xj,eflag)
        xi=z1
        f1=kappa_e2(iscont,e,xi,xj,eflag)	
      else
        xi=x1
        xj=z0
        f0=kappa_e2(iscont,e,xi,xj,eflag)
        xj=z1
        f1=kappa_e2(iscont,e,xi,xj,eflag)	
      endif
!
235   i=i+1
!
      if(i.gt.50)go to 237
        zm=z0-f0*(z1-z0)/(f1-f0)
!
        if(eixo.eq.1)then
          xj=x2
          xi=zm
          fm=kappa_e2(iscont,e,xi,xj,eflag)
        else
          xi=x1
          xj=zm
          fm=kappa_e2(iscont,e,xi,xj,eflag)
        endif	
!
        if(dabs(fm).lt.tol.or.dabs(z1-z0).lt.tol_in)then
          zk=zm
          go to 236
        else	
          gm=(f1-f0)/(z1-z0)
          hm=2.0d0*(f0+f1-2.0d0*fm)/((z1-z0)**2)
          del1=(-gm+dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
          del2=(-gm-dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
	
          if(dabs(del1).lt.dabs(del2))then
            del_min=del1
          else
            del_min=del2
          endif
!	
          zk=zm+del_min
	
          if(eixo.eq.1)then
              xi=zk
              xj=x2
              fk=kappa_e2(iscont,e,xi,xj,eflag)
          else
              xi=x1
              xj=zk
              fk=kappa_e2(iscont,e,xi,xj,eflag)
          endif	
!			
          if(dabs(fk).lt.tol.or.dabs(del_min).lt.tol_in)then
            go to 236
          else	
            if(del_min.gt.0.d0)then
	  z0=zm
	  f0=fm
	  z1=zk
	  f1=fk
            else
	  z0=zk
	  f0=fk
	  z1=zm
	  f1=fm
            endif
!
            goto 235
!
          endif	
      endif	

236   if(eixo.eq.1)then
        x1=zk
        x2=x2
      else
        x1=x1
        x2=zk
      endif	
!
237   return
      end

****************************************************************************
****************************************************************************
      SUBROUTINE bisec_ke(eixo,z_inf,z_sup,iscont,e,x1,x2,iter,eflag)
      double precision z_inf,z_sup,iscont,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj
      double precision kappa_e2,f0,f1,fm
      integer iter,eixo,maxiter,eflag
      external kappa_e2
      maxiter=100
!	
      if(e.lt.1.d-1)then
        tol=1.d-4
      else
        tol=1.d-5
      endif
c		
      iter=0
      z0=z_inf
      z1=z_sup
      if(eixo.eq.1)then
        xj=x2
        xi=z0
        f0=kappa_e2(iscont,e,xi,xj,eflag)
        xi=z1
        f1=kappa_e2(iscont,e,xi,xj,eflag)	
      else	
        xi=x1
        xj=z0
        f0=kappa_e2(iscont,e,xi,xj,eflag)
        xj=z1
        f1=kappa_e2(iscont,e,xi,xj,eflag)	
      endif
	
240   iter=iter+1	
      zm=(z0+z1)/2.0d0
c	
      if(eixo.eq.1)then
        xj=x2
        xi=zm
        fm=kappa_e2(iscont,e,xi,xj,eflag)
      else
        xi=x1
        xj=zm
        fm=kappa_e2(iscont,e,xi,xj,eflag)
      endif	
c	
      if(f0*fm.lt.0.d0.or.fm*f1.gt.0.d0)then
        z0=z0
        f0=f0
      else
        z0=zm
        f0=fm
      endif
c	
      if(dabs(fm).lt.tol)then
        go to 260	
      else
        if(iter.ge.maxiter)then
            zm=z0-f0*(z1-z0)/(f1-f0)
            go to 260
        else	
            go to 240
        endif
      endif
	
260   if(eixo.eq.1)then
        x1=zm
        x2=x2
      else
        x1=x1
        x2=zm
      endif
      RETURN
      END		
****************************************************************************
****************************************************************************
      SUBROUTINE raiz_lt(eixo,v_input,ks,e,x1,x2,eflag)
!       This subroutine find the root of the tangential eigenvalue of the 
!       mapping tensor lambda_t(ks,e,x1,x2)=0, keeping fixed x1 or x2 
!       Input Data:
!         eixo: eixo de calculo
!         	eixo=1, find the root along of the x1 axis
!        	eixo=2, find the root along of the x2 axis
!         	v_input: initial value 
!         	ks,e: parameters of the PNFW model
!       Output Data: the roots of this function
!
!       I/O Data
!         if eixo =1, x1(input),x2(output)
!         if eixo=2, x1(input),x2(output)				

      DOUBLE PRECISION v_input,ks,e,x1,x2,z_inf,z_sup
      INTEGER i,eixo,eflag
      call inter_lt(eixo,v_input,ks,e,x1,x2,z_inf,z_sup,eflag)
      call nrdf_lt(eixo,z_inf,z_sup,ks,e,x1,x2,i,eflag)
      if(i.gt.10)then
          write(*,*)'usando a bisecao'
          call bisec_lt(eixo,z_inf,z_sup,ks,e,x1,x2,i,eflag)	
      endif
!
      x1=dabs(x1)
      x2=dabs(x2)
!
      RETURN
      END
**********************************************************************************************
**********************************************************************************************
      subroutine inter_lt(eixo,v_input,ks,e,x1,x2,z_inf,z_sup,eflag)
      double precision v_input,ks,e,x1,x2,z_inf,z_sup
      double precision z,dz,xi,xj,lambda_t
      integer eixo,eflag
      EXTERNAL lambda_t
c	
      if(v_input.eq.0.d0)then
          z=1.d-3
      else
          z=v_input
      endif	
	
      dz=z*1.d-2
c	
      if(eixo.eq.1)then
          xi=z
          xj=x2
      else	
          xi=x1
          xj=z
      endif
c	
      if(lambda_t(ks,e,xi,xj,eflag).lt.0.d0)then
          z=z+dz
!
10        if(lambda_t(ks,e,xi,xj,eflag).lt.0.d0)then
	z_inf=z
	z=z+dz
	if(eixo.eq.1)then
	    xi=z
	    xj=x2
	else
	    xi=x1
	    xj=z
	endif				

	go to 10

          else
	z_sup=z
	if(eixo.eq.1)then
	    xi=z_sup
	    xj=x2
	else
	    xi=x1
	    xj=z
	endif				

	go to 30

          endif	

      else	
        z=z-dz
20        if(lambda_t(ks,e,xi,xj,eflag).gt.0.d0)then
	z_sup=z
	z=z-dz
	if(eixo.eq.1)then
	    xi=z
	    xj=x2
	else
	    xi=x1
	    xj=z
	endif				

	go to 20

          else
	z_inf=z
	if(eixo.eq.1)then
	    xi=z_inf
	    xj=x2
	else
	    xi=x1
	    xj=z_inf
	endif				

	go to 30
          endif			
      endif
c	
30    RETURN
      END
**********************************************************************************************
*********************************************************************************************
      subroutine nrdf_lt(eixo,z_inf,z_sup,ks,e,x1,x2,i,eflag)
      double precision z_inf,z_sup,ks,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj,zk,tol_in
      double precision lambda_t,f0,f1,fm,gm,hm
      double precision del1,del2,del_min,fk
      INTEGER i,eixo,eflag
      EXTERNAL lambda_t
c	
      i=0
      if(ks.lt.1.d-1)then
          tol=1.d-5
          tol_in=1.d-4
      else
          tol=1.d-6
          tol_in=1.d-5
      endif	
	
      z0=z_inf
      z1=z_sup
!
      if(eixo.eq.1)then
          xj=x2
          xi=z0
          f0=lambda_t(ks,e,xi,xj,eflag)
          xi=z1
          f1=lambda_t(ks,e,xi,xj,eflag)	
      else
          xi=x1
          xj=z0
          f0=lambda_t(ks,e,xi,xj,eflag)
          xj=z1
          f1=lambda_t(ks,e,xi,xj,eflag)	
      endif
				
35    i=i+1
      if(i.gt.50)go to 37
!
      zm=z0-f0*(z1-z0)/(f1-f0)
!	
      if(eixo.eq.1)then
          xj=x2
          xi=zm
          fm=lambda_t(ks,e,xi,xj,eflag)
      else
          xi=x1
          xj=zm
          fm=lambda_t(ks,e,xi,xj,eflag)
      endif	
!	
      if(dabs(fm).lt.tol.or.dabs(z1-z0).lt.tol_in)then
          zk=zm
          go to 36
      else	
!	
          gm=(f1-f0)/(z1-z0)
          hm=2.0d0*(f0+f1-2.0d0*fm)/((z1-z0)**2)
          del1=(-gm+dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
          del2=(-gm-dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
c	
          if(dabs(del1).lt.dabs(del2))then
              del_min=del1
          else
              del_min=del2
          endif
c	
          zk=zm+del_min
c	
          if(eixo.eq.1)then
              xi=zk
              xj=x2
              fk=lambda_t(ks,e,xi,xj,eflag)
          else
              xi=x1
              xj=zk
              fk=lambda_t(ks,e,xi,xj,eflag)
          endif	
c			
          if(dabs(fk).lt.tol.or.dabs(del_min).lt.tol_in)then
              go to 36
          else	
              if(del_min.gt.0.d0)then
	  z0=zm
	  f0=fm
	  z1=zk
	  f1=fk
              else
	  z0=zk
	  f0=fk
	  z1=zm
	  f1=fm
              endif
	  goto 35
          endif
      endif
		
36    if(eixo.eq.1)then
          x1=zk
          x2=x2
      else
          x1=x1
          x2=zk
      endif	
37    RETURN
      END
**********************************************************************************************
*********************************************************************************************	
      subroutine bisec_lt(eixo,z_inf,z_sup,ks,e,x1,x2,iter,eflag)
      double precision z_inf,z_sup,ks,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj
      double precision lambda_t,f0,f1,fm
      integer iter,eixo,maxiter,eflag
      external lambda_t
!
      maxiter=100
!	
      if(ks.lt.1.d-1)then
          tol=1.d-4
      else
          tol=1.d-5
      endif
!		
      iter=0

      z0=z_inf
      z1=z_sup
      if(eixo.eq.1)then
          xj=x2
          xi=z0
          f0=lambda_t(ks,e,xi,xj,eflag)
          xi=z1
          f1=lambda_t(ks,e,xi,xj,eflag)	
      else	
          xi=x1
          xj=z0
          f0=lambda_t(ks,e,xi,xj,eflag)
          xj=z1
          f1=lambda_t(ks,e,xi,xj,eflag)	
      endif
	
40    iter=iter+1	
      zm=(z0+z1)/2.0d0
!	
      if(eixo.eq.1)then
          xj=x2
          xi=zm
          fm=lambda_t(ks,e,xi,xj,eflag)
      else
          xi=x1
          xj=zm
          fm=lambda_t(ks,e,xi,xj,eflag)
      endif	
c	
      if(f0*fm.lt.0.d0.or.fm*f1.gt.0.d0)then
          z0=z0
          f0=f0
      else
          z0=zm
          f0=fm
      endif
c	
      if(dabs(fm).lt.tol)then
          go to 60	
      else
          if(iter.ge.maxiter)then
	zm=z0-f0*(z1-z0)/(f1-f0)
	go to 60
          else	
	go to 40
          endif
      endif
	
60    if(eixo.eq.1)then
          x1=zm
          x2=x2
      else
          x1=x1
          x2=zm
      endif
      RETURN
      END

****************************************************************************************
****************************************************************************************
      SUBROUTINE  raiz_lr(eixo,v_input,ks,e,x1,x2,eflag)
! This soubrutine find the root of the equation 
! lambda_r(ks,e,x1,x2)=0, keeping fixed x1 or  x2.
! Input Data:
! 	eixo: axis of calculation
! 	  eixo=1, calcula a raiz no eixo x1
! 	  eixo=2, calcula a raiz no eixo x2
! 	v_input: Input value for the calculation
! 	ks,e: PNFW lens parameters
! I/O Dats:
! 	se eixo =1, x1(saida),x2(entrada)
! 	se eixo=2, x1(entrada),x2(saida)				
!
      DOUBLE PRECISION v_input,ks,e,x1,x2,z_inf,z_sup
      INTEGER i,eixo,eflag
!
      call inter_lr(eixo,v_input,ks,e,x1,x2,z_inf,z_sup,eflag)	
      call nrdf_lr(eixo,z_inf,z_sup,ks,e,x1,x2,i,eflag)
      if(i.gt.10)then
          print*,'USING THE METHOD OF THE BISECTION'
          call bisec_lr(eixo,z_inf,z_sup,ks,e,x1,x2,i,eflag)
      endif	
!
      x1=dabs(x1)
      x2=dabs(x2)
!
      RETURN
      END
************************************************************************
************************************************************************
      SUBROUTINE inter_lr(eixo,v_input,ks,e,x1,x2,z_inf,z_sup,eflag)
      double precision v_input,ks,e,x1,x2,z_inf,z_sup
      double precision z,dz,xi,xj,lambda_r
      integer eixo,eflag
      EXTERNAL lambda_r
!	
      if(v_input.lt.1.d-5)then
          z=1.d-6
      else
          z=v_input
      endif	
!	
      dz=z*1.d-2
!	
      if(eixo.eq.1)then
          xi=z
          xj=x2
      else	
          xi=x1
          xj=z
      endif
!	
      if(lambda_r(ks,e,xi,xj,eflag).lt.0.d0)then
          z=z+dz
410       if(lambda_r(ks,e,xi,xj,eflag).lt.0.d0)then
!
              z_inf=z
              z=z+dz
!
              if(eixo.eq.1)then
	  xi=z
	  xj=x2
              else
	  xi=x1
	  xj=z
              endif				
!
              go to 410
!
          else
!
              z_sup=z
!
              if(eixo.eq.1)then
	  xi=z_sup
	  xj=x2
              else
	  xi=x1
	  xj=z
              endif				
!
              go to 430
          endif	
!
      else	
!
          z=z-dz
420       if(lambda_r(ks,e,xi,xj,eflag).gt.0.d0)then
              z_sup=z
              z=z-dz
              if(eixo.eq.1)then
	  xi=z
	  xj=x2
              else
	  xi=x1
	  xj=z
              endif				
!
              go to 420
!
          else
              z_inf=z
              if(eixo.eq.1)then
	  xi=z_inf
	  xj=x2
              else
	  xi=x1
	  xj=z_inf
              endif				
!
              go to 430
          endif			
      endif
!
430   RETURN
      END
*********************************************************
*********************************************************
      SUBROUTINE nrdf_lr(eixo,z_inf,z_sup,ks,e,x1,x2,i,eflag)
      double precision z_inf,z_sup,ks,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj,zk,tol_in
      double precision lambda_r,f0,f1,fm,gm,hm
      double precision del1,del2,del_min,fk
      INTEGER i,eixo,eflag
      EXTERNAL lambda_r
!	
      i=0
!
      if(ks.lt.1.d-1)then
          tol=1.d-4
          tol_in=1.d-4
      else
          tol=1.d-8
          tol_in=1.d-6
      endif
!	
      z0=z_inf
      z1=z_sup

      if(eixo.eq.1)then
          xj=x2
          xi=z0
          f0=lambda_r(ks,e,xi,xj,eflag)
          xi=z1
          f1=lambda_r(ks,e,xi,xj,eflag)	
      else
          xi=x1
          xj=z0
          f0=lambda_r(ks,e,xi,xj,eflag)
          xj=z1
          f1=lambda_r(ks,e,xi,xj,eflag)	
      endif
				
435   i=i+1
!
      if(i.gt.50)go to 437
!
      zm=z0-f0*(z1-z0)/(f1-f0)
      
      if(eixo.eq.1)then
          xj=x2
          xi=zm
          fm=lambda_r(ks,e,xi,xj,eflag)
      else
          xi=x1
          xj=zm
          fm=lambda_r(ks,e,xi,xj,eflag)
      endif	
!	
      if(dabs(fm).lt.tol.or.dabs(z1-z0).lt.tol_in)then
          zk=zm
          go to 436
!
      else	
!		Fazendo a aproximacao ate segunda ordem	
          gm=(f1-f0)/(z1-z0)
          hm=2.0d0*(f0+f1-2.0d0*fm)/((z1-z0)**2)
          del1=(-gm+dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
          del2=(-gm-dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
!	
          if(dabs(del1).lt.dabs(del2))then
              del_min=del1
          else
              del_min=del2
          endif
!	
          zk=zm+del_min
!	
          if(eixo.eq.1)then
              xi=zk
              xj=x2
              fk=lambda_r(ks,e,xi,xj,eflag)
          else
              xi=x1
              xj=zk
              fk=lambda_r(ks,e,xi,xj,eflag)
          endif	
c			
          if(dabs(fk).lt.tol.or.dabs(del_min).lt.tol_in)then
              go to 436
          else	
              if(del_min.gt.0.d0)then
	  z0=zm
	  f0=fm
	  z1=zk
	  f1=fk
              else
	  z0=zk
	  f0=fk
	  z1=zm
	  f1=fm
              endif

              goto 435

          endif
      endif
!
436   if(eixo.eq.1)then
          x1=zk
          x2=x2
      else
          x1=x1
          x2=zk
      endif
		
437   RETURN
      END
******************************************************************
******************************************************************	
      SUBROUTINE  bisec_lr(eixo,z_inf,z_sup,ks,e,x1,x2,iter,eflag)
      double precision z_inf,z_sup,ks,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj
      double precision lambda_r,f0,f1,fm
      integer iter,eixo,maxiter,eflag
      external lambda_r
!
      maxiter=100
!
	
      if(ks.lt.1.d-1)then
          tol=1.d-4
      else
          tol=1.d-5
      endif
c		
      iter=0
      z0=z_inf
      z1=z_sup
      if(eixo.eq.1)then
          xj=x2
          xi=z0
          f0=lambda_r(ks,e,xi,xj,eflag)
          xi=z1
          f1=lambda_r(ks,e,xi,xj,eflag)	
      else	
          xi=x1
          xj=z0
          f0=lambda_r(ks,e,xi,xj,eflag)
          xj=z1
          f1=lambda_r(ks,e,xi,xj,eflag)	
      endif
!	
440   iter=iter+1	
!
      zm=(z0+z1)/2.0d0
c	
      if(eixo.eq.1)then
          xj=x2
          xi=zm
          fm=lambda_r(ks,e,xi,xj,eflag)
      else
          xi=x1
          fm=lambda_r(ks,e,xi,xj,eflag)
      endif	
c	
      if(f0*fm.lt.0.d0.or.fm*f1.gt.0.d0)then
          z0=z0
          f0=f0
      else
          z0=zm
          f0=fm
      endif
c	
      if(dabs(fm).lt.tol)then
          go to 460	
      else
          if(iter.ge.maxiter)then
              zm=z0-f0*(z1-z0)/(f1-f0)
              go to 460
          else	
              go to 440
          endif
      endif
	
460   if(eixo.eq.1)then
          x1=zm
          x2=x2
      else
          x1=x1
          x2=zm
      endif
      RETURN
      END	
****************************************************************************
****************************************************************************
      SUBROUTINE  raiz_lw(eixo,lwu,v_input,ks,e,x1,x2,eflag)
! This soubroutine find the root of the Equation 
! L/W(ks,e,x1,x2)=(L/W)min, keeping fixed x1 or x2
! Input Data:
!        eixo: option of calculation
! 	eixo=1, calcula a raiz no eixo x1
! 	eixo=2, calcula a raiz no eixo x2
! 	lwt: Threshold value (L/W)_{min}
! 	If lwt<0, we have find L/W=-(L/W)_{min}
! 	Se lwt>0, we have find  L/W=(L/W)_{min}
! 	v_input: Input value 
! 	ks,e: parameters of the PNFW model
! I/ DATA:
!       if eixo =1, x1(input),x2(output)
!      if eixo= 2, x1(input),x2(output)				
!
      DOUBLE PRECISION lwu,lwt,v_input,ks,e,x1,x2,z_inf,z_sup
      INTEGER i,eixo,iter,eflag
!	
      lwt=1.0d0/lwu
!
      call inter_lw(eixo,lwt,v_input,ks,e,x1,x2,z_inf,z_sup,eflag)
      call nrdf_lw(eixo,lwt,z_inf,z_sup,ks,e,x1,x2,iter,eflag)
!
      if(iter.gt.10)then
          write(*,*)'USING THE METHOD OF THE BISECTION '
          call bisec_lw(eixo,lwt,z_inf,z_sup,ks,e,x1,x2,i,eflag)	
      endif	
!
      x1=dabs(x1)
      x2=dabs(x2)
!
      RETURN
      END
	
**********************************************************************************************
*********************************************************************************************
      subroutine inter_lw(eixo,lwt,v_input,ks,e,x1,x2,z_inf,z_sup,eflag)
      double precision lwt,v_input,ks,e,x1,x2,z_inf,z_sup
      double precision z,dz,xi,xj,wl
      integer eixo,eflag
      EXTERNAL wl
!	
      if(v_input.eq.0.d0)then
          z=1.d-3
      else
          z=v_input
      endif	
!	
      dz=z*1.d-3
	
      if(eixo.eq.1)then
          xi=z
          xj=x2
      else	
          xi=x1
          xj=z
      endif
c	
      if(wl(ks,e,xi,xj,eflag)-lwt.lt.0.d0)then
          z=z+dz
110       if(wl(ks,e,xi,xj,eflag)-lwt.lt.0.d0)then
              z_inf=z
              z=z+dz
              if(eixo.eq.1)then
	  xi=z
	  xj=x2
              else
	  xi=x1
	  xj=z
              endif
!
              go to 110
!
          else
!
              z_sup=z
              if(eixo.eq.1)then
	  xi=z_sup
	  xj=x2
              else
	  xi=x1
	  xj=z
              endif				
!
              go to 130
!
          endif	
      else	
!
          z=z-dz
120       if(wl(ks,e,xi,xj,eflag)-lwt.gt.0.d0)then
              z_sup=z
              z=z-dz
              if(eixo.eq.1)then
	  xi=z
	  xj=x2
              else
	  xi=x1
	  xj=z
              endif
!
              if(wl(ks,e,xi,xj,eflag)-lwt.lt.1.d-6)goto 121	
	  go to 120
              else
121	   z_inf=z
              if(eixo.eq.1)then
	  xi=z_inf
	  xj=x2
              else
	  xi=x1
	  xj=z_inf
              endif				
!
              go to 130
!
          endif			
      endif
!	
130   RETURN
      END	
**********************************************************************************************
*********************************************************************************************
      subroutine nrdf_lw(eixo,lwt,z_inf,z_sup,ks,e,x1,x2,i,eflag)
      double precision lwt,z_inf,z_sup,ks,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj,zk,tol_in
      double precision wl,f0,f1,fm,gm,hm
      double precision del1,del2,del_min,fk
      INTEGER i,eixo,eflag
      EXTERNAL wl
!	
      i=0
      if(ks.lt.1.d-1)then
          tol=1.d-4
          tol_in=1.d-5
      else
          tol=1.d-5
          tol_in=1.d-6
      endif
!		
      z0=z_inf
      z1=z_sup
      if(eixo.eq.1)then
          xj=x2
          xi=z0
          f0=wl(ks,e,xi,xj,eflag)-lwt
          xi=z1
          f1=wl(ks,e,xi,xj,eflag)-lwt	
      else
          xi=x1
          xj=z0
          f0=wl(ks,e,xi,xj,eflag)-lwt
          xj=z1
          f1=wl(ks,e,xi,xj,eflag)-lwt	
      endif		
!
135   i=i+1
!	
      if(i.gt.50) go to 137
!		
      zm=z0-f0*(z1-z0)/(f1-f0)
!	
      if(eixo.eq.1)then
          xj=x2
          xi=zm
          fm=wl(ks,e,xi,xj,eflag)-lwt
      else
          xi=x1
          xj=zm
          fm=wl(ks,e,xi,xj,eflag)-lwt
      endif
!		
      if(dabs(fm).lt.tol.or.dabs(z1-z0).lt.tol_in)then
          zk=zm
!
          go to 136
!
      else
          gm=(f1-f0)/(z1-z0)
          hm=2.0d0*(f0+f1-2.0d0*fm)/((z1-z0)**2)
          del1=(-gm+dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
          del2=(-gm-dsqrt(gm**2-4.0d0*hm*fm))/(2.0d0*hm)
!	
          if(dabs(del1).lt.dabs(del2))then
              del_min=del1
          else
              del_min=del2
          endif
c	
              zk=zm+del_min
c	
          if(eixo.eq.1)then
              xi=zk
              xj=x2
              fk=wl(ks,e,xi,xj,eflag)-lwt
          else
              xi=x1
              xj=zk
              fk=wl(ks,e,xi,xj,eflag)-lwt
          endif	
c			
          if(dabs(fk).lt.tol.or.dabs(del_min).lt.tol_in)then
              go to 136
          else	
              if(del_min.gt.0.d0)then
	  z0=zm
	  f0=fm
	  z1=zk
	  f1=fk
              else
	  z0=zk
	  f0=fk
	  z1=zm
	  f1=fk
              endif

	  goto 135

          endif	
      endif	
		
136   if(eixo.eq.1)then
          x1=zk
          x2=x2
      else
          x1=x1
          x2=zk
      endif	
137   RETURN
      END
**********************************************************************************************
**********************************************************************************************
      subroutine bisec_lw(eixo,lwt,z_inf,z_sup,ks,e,x1,x2,iter,eflag)
      double precision lwt,z_inf,z_sup,ks,e,x1,x2
      double precision tol,z0,z1,zm,xi,xj
      double precision wl,f0,f1,fm
      integer iter,eixo,maxiter,eflag
      external wl
!
      maxiter=100
      if(ks.lt.1.d-1)then
        tol=1.d-4
      else
        tol=1.d-5
      endif
c		
      iter=0
      z0=z_inf
      z1=z_sup
!
      if(eixo.eq.1)then
          xj=x2
          xi=z0
          f0=wl(ks,e,xi,xj,eflag)-lwt
          xi=z1
          f1=wl(ks,e,xi,xj,eflag)-lwt	
      else
          xi=x1
          xj=z0
          f0=wl(ks,e,xi,xj,eflag)-lwt
          xj=z1
          f1=wl(ks,e,xi,xj,eflag)-lwt	
      endif	
!	
140   iter=iter+1
! 	print*,'iteracao',iter
      zm=(z0+z1)/2.0d0
!
      if(eixo.eq.1)then
          xj=x2
          xi=zm
          fm=wl(ks,e,xi,xj,eflag)-lwt
      else
          xi=x1
          xj=zm
          fm=wl(ks,e,xi,xj,eflag)-lwt
      endif	
!	
      if(f0*fm.lt.0.d0.or.fm*f1.gt.0.d0)then
          z1=zm
          f1=fm
      else
          z0=zm
          f0=fm
      endif
	
      if(dabs(fm).lt.tol)then
          go to 160	
      else
          if(iter.gt.maxiter)then
              zm=z0-f0*(z1-z0)/(f1-f0)
              go to 160
          else	
              go to 140
          endif
      endif
	
160   if(eixo.eq.1)then
          x1=zm
          x2=x2
      else
          x1=x1
          x2=zm
      endif
c	
      RETURN
      END

****************************************************************************
****************************************************************************
