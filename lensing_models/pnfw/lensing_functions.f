! THIS FILE CONTAINS ALL THE BASIC LENSING FUNCTION OF THE PNFW LENS MODEL:
! ALMOST OF THESE FUNCTIONS ARE DEFINED IN Golse & Kneib A&A 390,821-827 (2002)
! FUNCTION f(x)     --> Eq. (4)
! FUNCTION g(x)     --> Eq. (5)
! FUNCTION alpha(x) --> Eq. 6(a)     It´s the deflection angle
! FUNCTION kappa(x)     --> Eq. 6(b)     It´s the characteristic convergence of the NFW model(*)
! FUNCTION gamma(x) --> Eq. 6(c)     It´s the shear of the NFW model(*)
! (*) These function are defined in the in the following functions
! FUNCTION kappa_e(x)     --> Eq.(17)     It´s the characteristic convergence of the PNFW model (**)
! FUNCTION gamma2_e(x)--> Eq.(18)     It´s the square shear of the PNFW model (**)
! (**) These function, here, aren´t independent of k_s
!
! AND ALSO SOME USEFUL FUNCTIONS
!
! FUNCTION lambda_r(x) : Radial eigenvalue of the mapping tensor
! FUNCTION lambda_t(x) : Tangential eigenvalue of the mapping thensor
! FUNCTION det_j(x): Is the absolute value of the determinat of the Mapping Tensor
! FUNCTION W/L : It´s the inverse of the lenght-to-width ratio  (L/W) of the eigenvalues
! FUNCTION kappa_e2 : This function is useful to calculate the iso-contours of the 
!	      convergence 
! FUNCTION y1_s: Mapping the ¨x2¨ points of the lens plane into the source plane points
! FUNCTION y2_s : Mapping the ¨x1¨ points  of the lens plane into the source plane points
!
      FUNCTION F(x)
!     Hint: arccosh(x)=ln(x+sqrt(x^2-1))
!     arcosh(1/x)=ln(1/x+sqrt(1/x^2-1))
!
      DOUBLE PRECISION x,xb,arg,f
!
      if(abs(x-1.).lt.1.d-6)then
        f=1.d0/3.d0
      endif
!	
      xb=x**2-1.d0
      if(x.lt.1.)then
        arg=1.d0/x+dsqrt(1.d0/(x**2)-1.d0)! argumento do arccosh(x)
        f=(1.d0/xb)*(1.d0-(1.d0/dsqrt(-xb))*dlog(arg))
      endif
!	
      if(x.gt.1.)then
        f=(1.d0/xb)*(1.d0-(1.d0/dsqrt(xb))*dacos(1.d0/x))
      end if
!
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION G(X)
!
      double precision d12,arg,x,g
!     Hint: arccosh(x)=ln(x+sqrt(x^2-1))
!	
      if(abs(x-1.d0).lt.1.d-5)then
        d12=1.d0/2.d0
        g=1.d0+dlog(d12)
      endif
!	
      if(x.lt.1.d0)then
        arg=1.d0/x+dsqrt(1.d0/(x**2)-1.d0)! argument of the arccosh(x)
        g=dlog(x/2.d0)+(1.d0/dsqrt(1.d0-x**2))*dlog(arg)
      endif
!	
      if(x.gt.1.d0)then
        g= dlog(x/2.d0) +(1.d0/dsqrt(x**2-1.d0))*dacos(1.d0/x)
      end if
!	
      RETURN
      END
****************************************************************************	
****************************************************************************      
      FUNCTION alpha(x)
!
      double precision alpha,x,g
      external g
      alpha=4.0d0*g(x)/x
!
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION kappa_e(e,x1,x2)
!     Input data
!     e : ellipticity of the angle deflection model
!     x1,x2: grid points

      double precision kappa_e,e,x1,x2
      double precision x1e,x2e,xe,dphi
      double precision f,g,kappa,gamma
      external f,g
!	
      x1e=x1*dsqrt(1.d0-e)
      x2e=x2*dsqrt(1.d0+e)
      xe=dsqrt(x1e**2+x2e**2)
!
      if(x1e.lt.dabs(1.d-8))then
        dphi=3.14159265d0
      else
        dphi=2.d0*datan(x2e/x1e)
      endif
!	
      kappa=2.d0*f(xe)
      gamma=2.d0*((2.d0*g(xe))/xe**2-f(xe))
!	
      kappa_e=kappa+e*dcos(dphi)*gamma
!	
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION gamma2_e(e,x1,x2)
!     Input Data 
!      e  ellipticity of the angle deflection method,
!      x1 e x2 grid points
!				
      double precision e,x1,x2,gamma2_e
      double precision x1e,x2e,xe,dphi
      double precision f,g,kappa,gamma
      external f,g
!	
      x1e=x1*dsqrt(1.d0-e)
      x2e=x2*dsqrt(1.d0+e)
      xe=dsqrt(x1e**2+x2e**2)
!
      if(x1e.lt.dabs(1.d-8))then
        dphi=3.14159265d0
      else
        dphi=2.d0*datan(x2e/x1e)
      endif
!	
      kappa=2.d0*f(xe)
      gamma=2.d0*((2.d0*g(xe))/xe**2-f(xe))
!	
      gamma2_e = gamma**2 +2.d0*e*dcos(dphi)*gamma*kappa +
     &           (kappa**2-(dsin(dphi)**2)*(gamma)**2)*e**2	
!     
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION lambda_t(ks,e,x1,x2)
!     Input Data
!     ks,e: parameters of the model
!     x1,x2: grid points		
!
      double precision lambda_t,ks,e,x1,x2
      double precision kappa_e,gamma2_e
      external kappa_e,gamma2_e
      lambda_t=1.d0-ks*(kappa_e(e,x1,x2)+dsqrt(gamma2_e(e,x1,x2)))
c	
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION lambda_r(ks,e,x1,x2)
!     Input Data
!     ks,e: parameters of the model
!     x1,x2: grid points		
!
      double precision lambda_r,ks,e,x1,x2
      double precision kappa_e,gamma2_e
      external kappa_e,gamma2_e
      lambda_r=1.d0-ks*(kappa_e(e,x1,x2)-dsqrt(gamma2_e(e,x1,x2)))
c	
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION det_j(ks,e,x1,x2)
!     Input Data
!     ks,e: parameters of the model
!     x1,x2: grid points	
      double precision det_j,ks,e,x1,x2
      double precision lambda_t,lambda_r
      external lambda_r,lambda_t
      det_j=dabs(lambda_t(ks,e,x1,x2)*lambda_r(ks,e,x1,x2))
!	
      RETURN	
      END	
****************************************************************************	
****************************************************************************
      FUNCTION  wl(ks,e,x1,x2)
!     W/L=lambda_t/lambda_r
!     Input Data:
!     ks,e: parameters of the model
!     x1,x2: grid points	
      double precision wl,ks,e,x1,x2
      double precision lambda_r, lambda_t
      external lambda_r,lambda_t
      wl=lambda_t(ks,e,x1,x2)/lambda_r(ks,e,x1,x2)
c	
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION kappa_e2(iscont,e,x1,x2)
!     Input Data
!     iscont: numerical value of the iso-contour
!     e : ellipticity of the angle deflection model
!     x1,x2: grid points
      double precision kappa_e2,e,x1,x2,iscont
      double precision kappa_e
      external kappa_e
c	
      kappa_e2=kappa_e(e,x1,x2)-iscont
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION	y1_s(ks,e,x1,x2)
!
      double precision y1_s,ks,e,x1,x2
      double precision x1e,x2e,xe,phie,alpha
      external alpha
      x1e=x1*dsqrt(1.d0-e)
      x2e=x2*dsqrt(1.d0+e)
      xe=dsqrt(x1e**2+x2e**2)
!
      if(x1e.lt.dabs(1.d-8))then
        phie=3.14159265d0/2.d0
      else
        phie= datan(x2e/x1e)
      endif
!
      y1_s=x1-ks*alpha(xe)*dsqrt(1.0d0-e)*dcos(phie)
!
      RETURN
      END
****************************************************************************	
****************************************************************************
      FUNCTION y2_s(ks,e,x1,x2)
!
      double precision y2_s,ks,e,x1,x2
      double precision x1e,x2e,xe,phie,alpha
      external alpha
      x1e=x1*dsqrt(1.d0-e)
      x2e=x2*dsqrt(1.d0+e)
      xe=dsqrt(x1e**2+x2e**2)
      if(x1e.lt.dabs(1.d-8))then
        phie=3.14159265d0/2.d0
      else
        phie=datan(x2e/x1e)
      endif
!
      y2_s=x2-ks*alpha(xe)*dsqrt(1.0d0+e)*dsin(phie)
!
      RETURN
      END
