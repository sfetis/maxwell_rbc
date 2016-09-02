 module get_aij

  contains
  
!!********************************************************************************************************
!! Процедура получения вектора правой части для заданного вектора облучения.
!!********************************************************************************************************
  
subroutine get_mid_RHSH(phi, jmax, rkt, norm1, norm2, tau, mod_k, f)  
  
    implicit none
    
    integer, intent(in) :: jmax
    real*8, intent(in) :: rkt(3, jmax), norm1(3, jmax), norm2(3, jmax),tau(3,2,jmax), mod_k, phi
    complex*16, intent(out) :: f (4 * jmax)
    
    integer :: i, j, m, l
    real*8 :: dot_kr, k(3)
    complex*16 :: e_ikr, vect_nxE0(3), E0(3)
      
      k(1) = cos(phi) * mod_k
      k(2) = sin(phi) * mod_k
      k(3) = 0d0
      
      m = 1
      
      do j = 1, jmax
      
	E0(1)= cmplx(sin(phi),0d0)
	E0(2)= cmplx(-cos(phi),0d0)
	E0(3)= cmplx(0d0,0d0)
	
	dot_kr = k(1) * rkt(1,j) + k(2) * rkt(2,j) + k(3) * rkt(3,j)
	
	e_ikr = EXP (cmplx(0d0,1d0) * dot_kr)
	
	E0(:) = E0(:) * e_ikr
	
	vect_nxE0 (1) =  (norm1(2,j) * E0(3) - norm1(3,j) * E0(2))
	vect_nxE0 (2) =  (norm1(3,j) * E0(1) - norm1(1,j) * E0(3))
	vect_nxE0 (3) =  (norm1(1,j) * E0(2) - norm1(2,j) * E0(1))	
		
	f(m) 		=(vect_nxE0(1) * tau(1,1,j)) + (vect_nxE0(2) * tau(2,1,j)) + (vect_nxE0(3) * tau(3,1,j))
	f(m + 1) 	=(vect_nxE0(1) * tau(1,2,j)) + (vect_nxE0(2) * tau(2,2,j)) + (vect_nxE0(3) * tau(3,2,j))  

	vect_nxE0 (1) =  (norm2(2,j) * E0(3) - norm2(3,j) * E0(2))
	vect_nxE0 (2) =  (norm2(3,j) * E0(1) - norm2(1,j) * E0(3))
	vect_nxE0 (3) =  (norm2(1,j) * E0(2) - norm2(2,j) * E0(1))
	
	f(2*jmax + m) 	=(vect_nxE0(1) * tau(1,1,j)) + (vect_nxE0(2) * tau(2,1,j)) + (vect_nxE0(3) * tau(3,1,j))
	f(2*jmax + m+1)	=(vect_nxE0(1) * tau(1,2,j)) + (vect_nxE0(2) * tau(2,2,j)) + (vect_nxE0(3) * tau(3,2,j))
	
	m = m + 2
	
      end do
      
  return
  end subroutine
  
subroutine get_mid_RHSV(phi, jmax, rkt, norm1, norm2, tau, mod_k, f)  
  
    implicit none
    
    integer, intent(in) :: jmax
    real*8, intent(in) :: rkt(3, jmax), norm1(3, jmax), norm2(3, jmax),tau(3,2,jmax), mod_k, phi
    complex*16, intent(out) :: f (4 * jmax)
    
    integer :: i, j, m, l
    real*8 :: dot_kr, k(3)
    complex*16 :: e_ikr, vect_nxE0(3), E0(3)
      
      k(1) = cos(phi) * mod_k
      k(2) = sin(phi) * mod_k
      k(3) = 0d0
      
      m = 1
      
      do j = 1, jmax
      
	E0(1)= cmplx(0d0,0d0)
	E0(2)= cmplx(0d0,0d0)
	E0(3)= cmplx(1d0,0d0)
	
	dot_kr = k(1) * rkt(1,j) + k(2) * rkt(2,j) + k(3) * rkt(3,j)
	
	e_ikr = EXP (cmplx(0d0,1d0) * dot_kr)
	
	E0(:) = E0(:) * e_ikr
	
	vect_nxE0 (1) =  (norm1(2,j) * E0(3) - norm1(3,j) * E0(2))
	vect_nxE0 (2) =  (norm1(3,j) * E0(1) - norm1(1,j) * E0(3))
	vect_nxE0 (3) =  (norm1(1,j) * E0(2) - norm1(2,j) * E0(1))	
		
	f(m) 		=(vect_nxE0(1) * tau(1,1,j)) + (vect_nxE0(2) * tau(2,1,j)) + (vect_nxE0(3) * tau(3,1,j))
	f(m + 1) 	=(vect_nxE0(1) * tau(1,2,j)) + (vect_nxE0(2) * tau(2,2,j)) + (vect_nxE0(3) * tau(3,2,j))  

	vect_nxE0 (1) =  (norm2(2,j) * E0(3) - norm2(3,j) * E0(2))
	vect_nxE0 (2) =  (norm2(3,j) * E0(1) - norm2(1,j) * E0(3))
	vect_nxE0 (3) =  (norm2(1,j) * E0(2) - norm2(2,j) * E0(1))
	
	f(2*jmax + m) 	=(vect_nxE0(1) * tau(1,1,j)) + (vect_nxE0(2) * tau(2,1,j)) + (vect_nxE0(3) * tau(3,1,j))
	f(2*jmax + m+1)	=(vect_nxE0(1) * tau(1,2,j)) + (vect_nxE0(2) * tau(2,2,j)) + (vect_nxE0(3) * tau(3,2,j))
	
	m = m + 2
	
      end do
      
  return
end subroutine
!!********************************************************************************************************
!! Процедура получения значения ЭПР в зависимости от вектора облучения
!!********************************************************************************************************

subroutine get_RCS(phi, jmax, JE, mod_k, tau, rkt, area, RCS)
  
    implicit none
    
    integer, intent(in) :: jmax
    
    complex*16, intent(in) :: JE(4 * jmax)
    
    real*8, intent(in) :: mod_k, tau(3,2,jmax), rkt(3, jmax), area(jmax)
    
    real*8, intent(out) :: RCS
    
    integer :: i, j, m
    
    complex*16 :: j_0(3), E0(3), e_ikr, dot_j_radar, tmp1(3), tmp(3), tmp2(3), j_m(3), i0
    
    real*8 :: exp2(3)
    
    real*8 :: phi, radar(3), dot_rad_x, dot_kr, pi, coef
    
      i0 = (0d0,1d0)
      coef = 1d0
      pi = 3.14159265d0
      radar(1) = cos(pi + phi)
      radar(2) = sin(pi + phi)
      radar(3) = 0d0

      tmp  = (0d0, 0d0)
      tmp1 = (0d0, 0d0)
      tmp2 = (0d0, 0d0)
      
      m = 1
      
      do j = 1, jmax
      
	j_0(:) = (JE(m) * tau(:,1,j) + JE(m+1) * tau(:,2,j))
	
	j_m(:) = (JE(m + 2*jmax) * tau(:,1,j) + JE((m+1)+2*jmax) * tau(:,2,j))
	
	dot_rad_x = radar(1) * rkt(1,j) + radar(2) * rkt(2,j) + radar(3) * rkt(3,j)
	
	dot_j_radar = j_0(1) * radar(1) + j_0(2) * radar(2) + j_0(3) * radar(3)
  
	e_ikr = exp(-(0d0,1d0) * mod_k * dot_rad_x)
	
	tmp2(1) = (radar(2) * j_m(3) - radar(3) * j_m(2))
	tmp2(2) = (radar(3) * j_m(1) - radar(1) * j_m(3))
	tmp2(3) = (radar(1) * j_m(2) - radar(2) * j_m(1))
	
	tmp2(:) = tmp2(:) * mod_k * i0
	
	tmp (:) = (j_0(:) - (radar(:)*dot_j_radar)) * mod_k * mod_k * i0 * coef
	
	tmp (:)  = (tmp(:) + tmp2(:)) * e_ikr  * area(j)	
	
	tmp1(:) = tmp1(:) + tmp(:)
	
	m = m + 2
	
	tmp(:) = 0d0
	tmp2(:)= 0d0
      enddo

	RCS = ((DIMAG(tmp1(1)))**2 + (DIMAG(tmp1(2)))**2 + (DIMAG(tmp1(3)))**2)
	RCS = RCS + ((REAL(tmp1(1)))**2 + (REAL(tmp1(2)))**2 + (REAL(tmp1(3)))**2)
	RCS = RCS * 4 * pi / 10000 
	RCS = 10 * log10(RCS)

    return
end subroutine
    
    
subroutine get_Ematr (rut0, rkt0, tau0, norm0, norm1, k, h_max, jmax, mtrx, m0, E_1, ker_1, ker_2)
  
  implicit none
    
    integer, intent(in) :: jmax, m0
    real*8, intent(in) :: rut0(3, 4, jmax), rkt0(3,jmax), norm0(3,jmax), tau0(3,2,jmax)
    real*8, intent(in) :: norm1(3, jmax), h_max, k
    complex*16, intent(out) :: mtrx(2*jmax, 2*jmax)
       
    external E_1, ker_1, ker_2
    
    integer :: i,j,p,q,m,l
    real*8 :: rb, coef
    complex*16 :: a(2,2), res(3), res2(3), res3(3)
    complex*16 :: i0

    
    rb	= 	2 * h_max/m0
    res2= 	(0d0,0d0)
    res3= 	(0d0,0d0)
    i0 = (0d0,1d0)
    coef = 1d0
    m = 1
    do i = 1, jmax
    
      l = 1
      res2 = 0d0
      res3 = 0d0
      
      do j = 1, jmax
	  
	  do p = 1, 2
	  
	    res2 = 0d0
	    res3 = 0d0
	    
	    call E_1(tau0(:,p,j), rkt0(:,i), rut0(:,:,j), rkt0(:,j), res2(:), rb, m0, m0, k, h_max, ker_2)

            call ker_1(rkt0(:,i), tau0(:,p,j), rut0(:,:,j), norm0(:,j), res3(:))   
            
	    res2(:) = res2(:) + res3(:)  
	    
            res(1) = (norm1(2,i) * res2(3) - norm1(3,i) * res2(2))
	    res(2) = (norm1(3,i) * res2(1) - norm1(1,i) * res2(3))
	    res(3) = (norm1(1,i) * res2(2) - norm1(2,i) * res2(1))
	    
	    do q = 1, 2
	    
	      a(q,p) = (res(1) * tau0(1,q,i) + res(2) * tau0(2,q,i) + res(3) * tau0(3,q,i))
	      
	    end do

	  end do
	  
      mtrx(m:(m+1),l:(l+1)) = a(:,:) *  coef * i0
      
      l = l + 2
      
      enddo
      
    m = m + 2
    
    enddo
  return
end subroutine get_Ematr
!   
subroutine get_Bmatr (rut0, rkt0, tau0, norm0, norm1, k, h_max, jmax, mtrx, m0, E_1, ker_1)
  
  implicit none
  
    integer, intent(in) :: jmax, m0
    real*8, intent(in) :: rut0(3, 4,jmax), rkt0(3,jmax), norm0(3,jmax), tau0(3,2,jmax)
    real*8, intent(in) :: norm1(3, jmax), h_max, k
    complex*16, intent(out) :: mtrx(2*jmax, 2*jmax)
       
    external E_1, ker_1
    
    integer :: i,j,p,q,m,l
    real*8 :: rb
    complex*16 :: a(2,2),res(3), res2(3), res3(3)
    complex*16 :: i0

    i0 = (0d0, 1d0)
    rb = 2 * h_max/m0
    res2 = (0d0,0d0)

    m = 1
    do i = 1, jmax
    
      l = 1
      res2 = 0d0
      
      do j = 1, jmax
	  
	  do p = 1, 2
	    res2 = 0d0
	    
	    call E_1(tau0(:,p,j), rkt0(:,i), rut0(:,:,j), rkt0(:,j), res2(:), rb, m0, m0, k, h_max, ker_1)
	    
            res(1) = (norm1(2,i) * res2(3) - norm1(3,i) * res2(2))
	    res(2) = (norm1(3,i) * res2(1) - norm1(1,i) * res2(3))
	    res(3) = (norm1(1,i) * res2(2) - norm1(2,i) * res2(1))
	    
	    do q = 1, 2
	    
	      a(q,p) = (res(1) * tau0(1,q,i) + res(2) * tau0(2,q,i) + res(3) * tau0(3,q,i))
	      
	    end do

	  end do
	  
      mtrx(m:(m+1),l:(l+1)) = - a(:,:)
      
      l = l + 2
      
      enddo
      
    m = m + 2
    
    enddo
  return
end subroutine get_Bmatr
  
subroutine get_db_matr(norm0, norm1, tau0, jmax, matr, mark)
  implicit none
  
  integer, intent(in)		:: jmax, mark
  complex*16, intent(inout)	:: matr(2*jmax, 2*jmax)
  real*8, intent(in)		:: norm1(3,jmax), norm0(3,jmax), tau0(3,2,jmax)
  integer			:: znak, l, k, i
  real*8 			:: pi
  complex*16 			:: res1(3), res2(3), res3(3), a(2,2)
  
  pi = 3.14159265d0
  
  if (mark==0) then
    znak = -1
  else if(mark==1) then
    znak = 1
  endif
  
  do i = 1, jmax
  
    do l = 1, 2
    
      res1(1) = norm0(2,i) * tau0(3,l,i) - norm0(3,i) * tau0(2,l,i)
      res1(2) = norm0(3,i) * tau0(1,l,i) - norm0(1,i) * tau0(3,l,i)
      res1(3) = norm0(1,i) * tau0(2,l,i) - norm0(2,i) * tau0(1,l,i)  
      
      res2(1) = norm1(2,i) * res1(3) - norm1(3,i) * res1(1)
      res2(2) = norm1(3,i) * res1(1) - norm1(1,i) * res1(3)
      res2(3) = norm1(1,i) * res1(2) - norm1(2,i) * res1(1)   
      
      do k = 1, 2
	a(k,l) = res2(1) * tau0(1,k,i) + res2(2) * tau0(2,k,i) + res2(3) * tau0(3,k,i)
      enddo
  
    enddo
      matr((2*i-1):2*i,(2*i-1):2*i) = matr((2*i-1):2*i,(2*i-1):2*i) + a(:,:)*znak * 2 * pi 
  enddo
end subroutine get_db_matr
    
subroutine get_da_matr( norm0, norm1, rut0, tau0, area0, jmax, matr, k0, soot, ind, mark)
  
  integer, intent(in)		:: jmax, mark, soot(4, jmax), ind(jmax)
  complex*16, intent(inout)	:: matr(2*jmax, 2*jmax)
  real*8, intent(in)		:: norm1(3,jmax), norm0(3,jmax), tau0(3,2,jmax), rut0(3,4,jmax), area0(jmax), k0
  integer			:: znak, l, k, m, i, j
  real*8 			:: r(3,4), dot_j_rn, pi, coef
  complex*16 			:: a(2,2),res1(3), res2(3), res3(3), i0
  pi = 3.14159265d0
  i0 = (0d0,1d0)
  coef = 1d0
  if (mark == 0) then
    znak = 1

  else if(mark==1) then
    znak = -1

  endif
  
  do i = 1, jmax
  
  call get_Rm(rut0(:,:,i), r(:,:))
  
    do m = 1, ind(i)
    
      j = soot(m, i)
      
      if (j/=0) then

      res1(1) = norm1(2,i) * norm0(3,i) - norm1(3,i) * norm0(2,i)
      res1(2) = norm1(3,i) * norm0(1,i) - norm1(1,i) * norm0(3,i)
      res1(3) = norm1(1,i) * norm0(2,i) - norm1(2,i) * norm0(1,i)
      
      res2(1) = r(2,m) * norm0(3,i) - r(3,m) * norm0(2,i)
      res2(2) = r(3,m) * norm0(1,i) - r(1,m) * norm0(3,i)
      res2(3) = r(1,m) * norm0(2,i) - r(2,m) * norm0(1,i)
      
      do k = 1, 2
	do l = 1, 2
	  dot_j_rn = res2(1) * tau0(1,l,j) + res2(2) * tau0(2,l,j) + res2(3) * tau0(3,l,j)
	  
	  res1(:) = res1(:) * dot_j_rn/2d0
	  
	  a(k,l) = res1(1) * tau0(1,k,i) + res1(2) * tau0(2,k,i) + res1(3) * tau0(3,k,i)	  
	enddo
      enddo
      
      matr((2*i-1):2*i, (2*j-1):2*j) = matr((2*i-1):2*i,(2*j-1):2*j) + a(:,:) * znak * 2d0 * pi * i0
      
      end if  
    enddo
  enddo
end subroutine
    
subroutine get_Rm(rut,r)
  implicit none
  
  real*8, intent(in) :: rut(3,4)
  real*8, intent(out) :: r(3,4)
  integer :: i
  
  do i = 1, 3
    r(:,i) = rut(:,4-i) - rut(:,4-i+1)
  enddo 
    r(:,4) = rut(:,4) - rut(:,1)
end subroutine   
  
subroutine get_Ex(jmax, je, jm, x, rut, rkt, norm, area, k, hmax, m0, res, soot, ind, znak, E_1, ker1, ker2, ker3, rm)
    implicit none
    
    integer, intent(in)		:: jmax, soot(4, jmax), ind(jmax), znak, m0
    complex*16, intent(in)	:: je(3, jmax), jm(3, jmax)
    real*8, intent(in) 		:: rut(3, 4, jmax), rkt(3, jmax), norm(3, jmax), x(3), area(jmax), k, hmax
    external ker1, ker2, ker3, E_1, rm
    complex*16, intent(out)	:: res(3)
    
    complex*16			:: res1(3), res2(3), res3(3), res4(3), res5(3), tmp(3), tmp2(3), tmp1(3), i0
    integer 			:: i, m, l, j
    real*8			:: r(3,4), pi
    
    pi = 3.1415926d0
    i0 = (0d0, 1d0)

    do i = 1, jmax
      call E_1(je(:,i), x, rut(:,:,i), rkt(:,i), res1, 2*hmax/m0, m0, m0, k, hmax, ker1)
      
      call ker2(x, je(:,i), rut(:,:,i), norm(:,i),res2)
      
      res1 = (res1 + res2) * i0
      
      call E_1(jm(:,i), x, rut(:,:,i), rkt(:,i), res3, 2*hmax/m0, m0, m0, k, hmax, ker3)
      
      
      
      res4(1) = (norm(2,i) * jm(3,i) - norm(3,i) * jm(2,i))
      res4(2) = (norm(3,i) * jm(1,i) - norm(1,i) * jm(3,i))
      res4(3) = (norm(1,i) * jm(2,i) - norm(2,i) * jm(1,i))
      
      res4(:) = - znak * res4(:) * 2 * pi
      
      call rm(rut(:,:,i), r(:,:))
      res5(:) = (0d0,0d0)
      tmp(:) = (0d0,0d0)
      do l = 1, ind(i)
	
	j = soot(l,i)
	if (j/=0) then
	
	tmp1(1) = r(2,l)*norm(3,i) - r(3,l)*norm(2,i)
	tmp1(2) = r(3,l)*norm(1,i) - r(1,l)*norm(3,i)
	tmp1(3) = r(1,l)*norm(2,i) - r(2,l)*norm(1,i)
	
	tmp(:) = (je(:, j) + je(: ,i))/2d0
	
	res5 = res5 + tmp1(1)*tmp(1) + tmp1(2)*tmp(2) + tmp1(3)*tmp(3)
	
	endif
      enddo
      
      res5 = znak * res5 * 2 * pi * i0  / area(i)
      
      res = res1 + res3 + res4 + res5 + res
    enddo
end subroutine get_Ex
end