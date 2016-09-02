module AV_sub
  contains
!***************************************************************
!   
!***************************************************************
	subroutine E_1(g, x, rut0, rkt0, e, rb, m_0, n_0, k, h_max, K_0)
!*************************************************************

implicit none
  real*8 x(3), rut0(3,4), rb, k, h_max, rkt0(3), g(3)
  integer j,i, m, n, m_0, n_0
  complex*16 e(3), e1(3)
external K_0

real*8 p,q,p1,q1, a(3), b(3), a1(3),a2(3),a3(3), a4(3),m1(3),m2(3),rn(3),s,rc(3)

real*8 dr0,q0


dr0 = sqrt( (x(1)-rkt0(1))**2 + (x(2)-rkt0(2))**2 + (x(3)-rkt0(3))**2 )

q0 = dr0/h_max

m=1; n=1

if( q0.le.10) then; m=m_0/4; n=n_0/4; end if
if( q0.le.6) then; m=m_0/2; n=n_0/2; end if
if( q0.le.3) then; m=m_0; n=n_0; end if
if(m.lt.1) m=1; if(n.lt.1) n=1; 

!m=20; n=20

e = 0

  do i=1,m
  do j=1,n
		p = (i-1.)/m
		q = (j-1.)/n
		p1 = (i-0.)/m
		q1 = (j-0.)/n

		a = q*rut0(:,2) + (1-q)*rut0(:,1); b = q*rut0(:,3) + (1-q)*rut0(:,4)
		a1 = p*b + (1-p)*a
		a4 = p1*b + (1-p1)*a
		a = q1*rut0(:,2) + (1-q1)*rut0(:,1); b = q1*rut0(:,3) + (1-q1)*rut0(:,4); 
		a2 = p*b + (1-p)*a
		a3 = p1*b + (1-p1)*a

		rc = (a1 + a2 + a3 + a4) / 4.

		m1 = ((a2 + a3) - (a1 + a4)) / 2.
		m2 = ((a3 + a4) - (a1 + a2)) / 2.
		rn(1) = m1(2)*m2(3) - m1(3)*m2(2) 
		rn(2) = m1(3)*m2(1) - m1(1)*m2(3)
		rn(3) = m1(1)*m2(2) - m1(2)*m2(1)

		s = sqrt(rn(1)**2 + rn(2)**2 + rn(3)**2)
		
		call K_0(g, x, rc, k, e1, rb)
		

		e = e + e1 * s

  end do
  end do		

return
end subroutine

subroutine E_2(g, x, rut0, rkt0, e, rb, m_0, n_0, k, h_max, K_0)
!*************************************************************

implicit none
  real*8 x(3), rut0(3,4), rb, k, h_max, rkt0(3)
  integer j,i, m, n, m_0, n_0
  complex*16 e(3), e1(3), g(3)
external K_0

real*8 p,q,p1,q1, a(3), b(3), a1(3),a2(3),a3(3), a4(3),m1(3),m2(3),rn(3),s,rc(3)

real*8 dr0,q0


dr0 = sqrt( (x(1)-rkt0(1))**2 + (x(2)-rkt0(2))**2 + (x(3)-rkt0(3))**2 )

q0 = dr0/h_max

m=1; n=1

if( q0.le.10) then; m=m_0/4; n=n_0/4; end if
if( q0.le.6) then; m=m_0/2; n=n_0/2; end if
if( q0.le.3) then; m=m_0; n=n_0; end if
if(m.lt.1) m=1; if(n.lt.1) n=1; 

!m=20; n=20

e = 0

  do i=1,m
  do j=1,n
		p = (i-1.)/m
		q = (j-1.)/n
		p1 = (i-0.)/m
		q1 = (j-0.)/n

		a = q*rut0(:,2) + (1-q)*rut0(:,1); b = q*rut0(:,3) + (1-q)*rut0(:,4)
		a1 = p*b + (1-p)*a
		a4 = p1*b + (1-p1)*a
		a = q1*rut0(:,2) + (1-q1)*rut0(:,1); b = q1*rut0(:,3) + (1-q1)*rut0(:,4); 
		a2 = p*b + (1-p)*a
		a3 = p1*b + (1-p1)*a

		rc = (a1 + a2 + a3 + a4) / 4.

		m1 = ((a2 + a3) - (a1 + a4)) / 2.
		m2 = ((a3 + a4) - (a1 + a2)) / 2.
		rn(1) = m1(2)*m2(3) - m1(3)*m2(2) 
		rn(2) = m1(3)*m2(1) - m1(1)*m2(3)
		rn(3) = m1(1)*m2(2) - m1(2)*m2(1)

		s = sqrt(rn(1)**2 + rn(2)**2 + rn(3)**2)
		
		call K_0(g, x, rc, k, e1, rb)
		

		e = e + e1 * s

  end do
  end do		

return
end subroutine
end
 