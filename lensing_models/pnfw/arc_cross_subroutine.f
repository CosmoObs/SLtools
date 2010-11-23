      subroutine c_section(nptx,v_input,ks,e,lw_u,muth,iflag,v_out,
     &	sigma_q,sigma_u,u_min)
c     This subroutine calculate:
c     If  iflag=2  	
c     	-Dimensionless arc cross section
c     If iflag=3
c	-Dimensionless arc cross section with u>u_th constrain
c     Variables definition
c	iv_lt: Input valu for the lambda_t root
c	x1_lt: intersection of the lambda_t curve with the axis x1
c               And similarly to ql=-10 (iv_lwn,...), q_l=10 (iv_lwp)
c               
c	sum_x2: Area in the lens plane (axis x2) (first quadrant)
c	s_area : Area in the source plane limited by the curves  |q_l|>=q_l,min (first quadrant)
c	sigma_q: Total area in the source plane 
c               sigma_u: Total area in the source lens taking into account u>u_th
c	External subroutines:
c	subroutine raiz_**(method,eixo,v_input,(lw_u)ks,e,x1,x2)
c	** corresponding to lambda_t ou  L/W ratio
c	method (1) Newton-Rapshon Method 
c	       (2) Bisection Method 
c	eixo (1) Calculating the root in the x1 axis for some value of  x2
c	     (2) Calculating the root in the x2 axis for some values of  x1 
c	I/O variables
      double precision v_input,v_in,ks,e,lw_u,muth
	double precision v_out,sigma_q,sigma_u,u_min
c	Points and steps of the adaptative grid 
      double precision x1,x2,dx1,dx2,x2_u
c	Variables related with lambda_t
      double precision iv_lt,x1_lt
c	Variables related with  q_l=-10
      double precision iv_lwn,x1q_n,x2q_n
c	Variables related with  q_l=10			
      double precision iv_lwp,x1q_p,x2q_p
c	External functions
      double precision det_j,wl,s_area,sum_x2
      double precision det_j2,det_jmax
      double precision x1m,x2p,x2n
      integer nptx,iflag,i,j,eixo,npt1,npt2
      dimension x1m(nptx),x2p(nptx),x2n(nptx)
      external det_j,wl
c     Defining the input values
      
C
       v_in=v_input
c
c     Calculating the intersections of the R_lambda curves
      iv_lt=v_in
      eixo=1
      x2=0.d0
c     Calculating the intersection R_lambda=infty curve with x1 axis
      call raiz_lt(eixo,iv_lt,ks,e,x1_lt,x2)
cc    Calculating the intersection of the R_l> 0 with x1 axis
      iv_lwp=x1_lt
      call raiz_lw(eixo,lw_u,iv_lwp,ks,e,x1q_p,x2)
      x1m(nptx)=x1q_p
      x2p(nptx)=0.0d0
      x2n(nptx)=0.d0
c     Calculating the intersection of   R_l<0	with x1 axis
      iv_lwn=x1_lt
      call raiz_lw(eixo,-lw_u,iv_lwn,ks,e,x1q_n,x2)
      write(*,*)'intersections :',x1q_n,x1_lt,x1q_p
      npt1=nptx
      dx1=(x1q_p)/dfloat(npt1-1)
c
      v_out=x1_lt
      s_area=0.0d0
      iv_lwp=0.d0
      iv_lwn=0.d0
      do i=npt1-1,1,-1
          x1=0.d0+dx1*(i-1)
          x1m(i)=x1
          call raiz_lw(2,lw_u,iv_lwp,ks,e,x1,x2q_p)
          x2p(i)=x2q_p
          iv_lwp=x2q_p
          if(x1.gt.x1q_n)then
              x2q_n=0.d0
              x2n(i)=x2q_n
          else
              iv_lt=x2q_p
              call raiz_lt(2,iv_lt,ks,e,x1,x2)
              iv_lwn=x2	
              call raiz_lw(2,-lw_u,iv_lwn,ks,e,x1,x2q_n)
              x2n(i)=x2q_n
              iv_lwn=x2q_n
          endif   
c     Calculating the sum along the x2 axis
          if(x1.ge.x1q_n)then
              npt2=nptx
          else
              npt2=(nptx+1)/2
          endif  
          dx2=(x2q_p-x2q_n)/dfloat(npt2-1)
          sum_x2=0.0d0
          do j=1,npt2
              x2=x2q_n+(j-1)*dx2
              sum_x2=sum_x2+det_j(ks,e,x1,x2)*dx2
          end do
          if(i.gt.1)s_area=s_area+dx1*sum_x2
      end do    
      write(*,*)'area in first quadrant of the source plane ',s_area
      sigma_q=4.0d0*s_area
      write(*,*)'area in the source plane ',sigma_q
c     Finding the minimum value of the magnification
      npt1=nptx	
      det_j2=0.d0
      do i=1,npt1
          det_j1=det_j(ks,e,x1m(i),x2p(i))
          if(det_j1.gt.det_j2)then
              det_j2=det_j1
          endif
      end do
      det_jmax=det_j2	
      u_min=1.0d0/det_jmax 
	write(*,*)'minimum magnification ',u_min
      sigma_u=0.d0
      if(iflag.eq.2)goto 100
c     Calculating the dimensionless arc cross section with the constraint u>u_th
      npt1=nptx
      dx1=(x1q_p)/dfloat(npt1-1)
      s_area=0.d0
      do i=2,npt1
          if(x1m(i).lt.x1q_n)then
              npt2=nptx
          else
              npt2=(nptx+1)/2
          endif
          dx2=(x2p(i)-x2n(i))/dfloat(npt2-1)
          sum_x2=0.d0
          do j=1,npt2
              x2_u=x2n(i)+(j-1)*dx2
	if(det_j(ks,e,x1m(i),x2_u).le.1.0d0/muth)then !(si mu >= mu_th)
	    sum_x2=sum_x2+det_j(ks,e,x1m(i),x2_u)*dx2
	endif
          end do
          s_area=s_area+sum_x2*dx1
      end do
      sigma_u=4.0d0*s_area	
100   return
      end
	