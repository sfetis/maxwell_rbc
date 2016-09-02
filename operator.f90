module operators_int
  
  implicit none
  contains

!!*************************************************************************
!! Процедура, которая вычисляет значение оператора K_0
!!*************************************************************************
  subroutine kernel_00 (x, g, cell, norm, res)
    implicit none
    
    real*8, intent (in) :: g(3)
    complex*16, intent (out):: res(3)
    real*8, intent (in) :: x(3), cell(3,4), norm(3)
    complex*16 :: temp(3)
    integer :: i
    
    res = (0d0, 0d0)
    
    call lin_int_Km0(x, cell(:,4), cell(:, 1), norm, g, res)
    
    do i = 1, 3
    
      call lin_int_Km0(x, cell(:, i), cell(:, i + 1), norm, g, temp)
      
      res = res + temp
      
    enddo
    return
  end subroutine kernel_00
  
!!************************************************************************* 
!! Вычисление К_0 по отрезку
!!************************************************************************* 
 
  subroutine lin_int_Km0(x, a, b, norm, g, res)
    implicit none
    real*8, intent(in) 	:: g(3)
    complex*16, intent(out) 	:: res(3)
    real*8, intent(in)		:: x(3), a(3), b(3), norm(3)
    real*8, dimension(3)	:: x_a, x_b, b_a, tau
    real*8			:: mod_xa, mod_xb, mod_ba, dot_xa_xb
    complex*16 			:: vect_jn(3),  dot_jn_tau
    
    x_a = x - a
    x_b = x - b
    b_a = a - b
    
    mod_xa = sqrt(x_a(1)**2 + x_a(2)**2 + x_a(3)**2)
    mod_xb = sqrt(x_b(1)**2 + x_b(2)**2 + x_b(3)**2)
    mod_ba = sqrt(b_a(1)**2 + b_a(2)**2 + b_a(3)**2)
    
    dot_xa_xb = x_a(1) * x_b(1) + x_a(2) * x_b(2) + x_a(3) * x_b(3)
    
    vect_jn(1) = g(2) * norm(3) - g(3) * norm(2)
    vect_jn(2) = g(3) * norm(1) - g(1) * norm(3)
    vect_jn(3) = g(1) * norm(2) - g(2) * norm(1)
    
    tau =  b_a
    
    dot_jn_tau = vect_jn(1) * tau(1) + vect_jn(2) * tau(2) + vect_jn(3) * tau(3)
    
    res = dot_jn_tau * ( (x_b/mod_xb) + (x_a/mod_xa) ) / ( (mod_xa*mod_xb) + dot_xa_xb)

    return
  end subroutine lin_int_Km0

  
!!*************************************************************************  
!! Процедура, которая вычисляет подынтегральную функцию для оператора K_1
!! Формула взята из отчёта за 2013 год.
!!*************************************************************************  
  
subroutine kernel_2 (g, x, y, k, res, eps)
    implicit none
  
    real*8, intent(in) :: x(3), y(3), k, eps
    real*8, intent(in) :: g(3)
    complex*16, intent(out) :: res(3)
    complex*16 :: i, dot_rg, e_ikr
    real*8 :: th, r(3), mod_r
    logical :: flag_1
    
    th = 1d0
    i = (0d0, 1d0)
    r =  x - y
    mod_r = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    
    r = r / mod_r
    
    dot_rg = (r(1)*g(1) + r(2)*g(2) + r(3)*g(3))
    
    e_ikr = exp(i * k * mod_r)
    
    flag_1 = .false.

    if (mod_r < 1d-8) then
      flag_1 = .true.
    else
      if (mod_r < eps) then
	th = (3 * ((mod_r/eps)**2)) - (2 * (mod_r/eps)**3) 
      else
	th = 1d0	
      end if
      
    end if
    
    if (flag_1 .eqv. .true.) then
      res = 0d0
    else
    
      res = g * ((i*k/mod_r**2) + ((k**2)/mod_r))
      res = res + r * dot_rg * (((- 3 * i * k)/(mod_r**2)) - (k**2/mod_r))
      res = res * exp(i * k * mod_r)
      res = res + (((1 - exp(i*k*mod_r))/mod_r**3)*(g - 3*r*dot_rg))
      res = res * th
      
    endif
    return
end subroutine kernel_2

subroutine kernel_1(g, x, y, k, res, eps)
    implicit none
  
    real*8, intent(in) :: x(3), y(3), k, eps
    complex*16, intent(in) :: g(3)
    complex*16, intent(out) :: res(3)
    complex*16 :: i, dot_rg, e_ikr
    real*8 :: th, r(3), mod_r
    logical :: flag_1
    
    th = 1d0
    i = (0d0, 1d0)
    r =  x - y
    mod_r = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    r = r / mod_r
    
    dot_rg = (r(1)*g(1) + r(2)*g(2) + r(3)*g(3))
    
    e_ikr = exp(i * k * mod_r)
    
    flag_1 = .false.

    if (mod_r < 1d-10) then
      flag_1 = .true.
    else
    
      if (mod_r < eps) then
	th = (3 * ((mod_r/eps)**2)) - (2 * (mod_r/eps)**3) 
      else
	th = 1d0
      end if
      
    end if
    
    if (flag_1 .eqv. .true.) then
      
      res = 0d0
    else
      res = (g - 3 * r * dot_rg) * ( 1 - exp( i * k * mod_r)  + i * k * mod_r * exp( i * k * mod_r)) / (mod_r**3)
      res = res + ((g - r * dot_rg) * k * k * exp( i * k * mod_r) / mod_r)
      res = res * th
    endif
    return
  end subroutine kernel_1  
  
  subroutine kernel_3(g, x, y, k, res, eps)
    implicit none
  
    real*8, intent(in)		::k, eps
    real*8, intent(in)		::g(3), x(3), y(3)
    complex*16, intent(out)	::res(3)
    complex*16			::i, res2(3)
    real*8 			::th, r(3), mod_r
    logical 			::flag_1
    
    i = (0d0, 1d0)
    r = y - x
    mod_r = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    
    flag_1 = .false.

    if (mod_r < 1d-10) then
      flag_1 = .true.
    else
    
      if (mod_r < eps) then
	th = (3 * (mod_r/eps)**2) - (3 * (mod_r/eps)**2) 
      else
	th = 1d0
      end if
      
    end if
    
    if (flag_1 .eqv. .true.) then
      res = 0d0
    else

      res2 = (1d0, 0d0) * r
      res2 = res2 * (i * k * mod_r - 1d0) * th
      res2 = res2 * exp( i * k * mod_r) / (mod_r**3)

      res(1) = (g(2) * res2(3) - g(3) * res2(2))
      res(2) = (g(3) * res2(1) - g(1) * res2(3))
      res(3) = (g(1) * res2(2) - g(2) * res2(1))
      
    end if
    return
  end subroutine kernel_3
  
  
  subroutine get_soot(jmax, h, l, soot, ind)
    implicit none
    integer, intent(in) :: jmax, l, h
    integer, intent(out) :: soot(4,jmax), ind(jmax)
    integer j
    
    soot = 0
    ind = 0
    
! угловые точки - два соседа

    ind(1) = 2
    ind(l) = 2
    ind(jmax) = 2
    ind (jmax-l+1) = 2
    
    soot(1,1) = 2
    soot(2,1) = 1 + l
    
    soot(2,l) = 2 * l
    soot(3,l) = l - 1
    
    soot(3,jmax) = jmax - 1
    soot(4,jmax) = jmax - l
    
    soot(4,jmax - l + 1) = jmax - 2*l + 1
    soot(1,jmax - l + 1) = jmax - l + 2
    
! боковые стенки вдоль размаха - три соседа

    do j = 2, l-1
    
      ind(j) = 3
      
      soot(1, j) = j + 1
      soot(2, j) = j + l
      soot(3, j) = j - 1
      
      ind(jmax - l + j) = 3
				  
      soot(1, jmax - l + j ) = jmax - l + j - 1
      soot(3, jmax - l + j ) = jmax - l + j + 1
      soot(4, jmax - l + j ) = jmax - l + j - l
    enddo

! боковые стенки вдоль хорды - три соседа
    
    do j = 1, h-2
    
      ind (j * l + 1) = 3
      soot(1, j * l + 1) = j * l + 2
      soot(2, j * l + 1) = j * l + 1 + l
      soot(4, j * l + 1) = j * l + 1 - l
      
    enddo 
    
    do j = 2, h-1
    
      ind (j * l) = 3
      soot(2, j * l) = j * l + l
      soot(3, j * l) = j * l - 1
      soot(4, j * l) = j * l - l
      
    enddo
    
    do j = 1, jmax
    
      if(ind(j)==0) then
	ind(j) = 4
	soot(1, j) = j + 1 
	soot(2, j) = j + l 
	soot(3, j) = j - 1
	soot(4, j) = j - l
      endif
      
    enddo
    
    end subroutine get_soot 

    
   subroutine get_soot_el(jmax, h, l, soot, ind)
    implicit none
    integer, intent(in) :: jmax, l, h
    integer, intent(out) :: soot(4,jmax), ind(jmax)
    integer j
    
    soot = 0
    ind = 0
    
! угловые точки - два соседа

    ind(1) = 2
    ind(l) = 2
    ind(jmax) = 2
    ind (jmax-l+1) = 2
    
    soot(1,1) = 2
    soot(2,1) = 1 + l
    
    soot(2,l) = 2 * l
    soot(3,l) = l - 1
    
    soot(3,jmax) = jmax - 1
    soot(4,jmax) = jmax - l
    
    soot(4,jmax - l + 1) = jmax - 2*l + 1
    soot(1,jmax - l + 1) = jmax - l + 2
    
! боковые стенки вдоль размаха - три соседа

    do j = 2, l-1
    
      ind(j) = 3
      
      soot(1, j) = j + 1
      soot(2, j) = j + l
      soot(3, j) = j - 1
      
      ind(jmax - l + j) = 3
				  
      soot(1, jmax - l + j ) = jmax - l + j - 1
      soot(3, jmax - l + j ) = jmax - l + j + 1
      soot(4, jmax - l + j ) = jmax - l + j - l
    enddo

! боковые стенки вдоль хорды - три соседа
    
    do j = 1, h-2
    
      ind (j * l + 1) = 3
      soot(1, j * l + 1) = j * l + 2
      soot(2, j * l + 1) = j * l + 1 + l
      soot(4, j * l + 1) = j * l + 1 - l
      
    enddo 
    
    do j = 2, h-1
    
      ind (j * l) = 3
      soot(2, j * l) = j * l + l
      soot(3, j * l) = j * l - 1
      soot(4, j * l) = j * l - l
      
    enddo
    
    do j = 1, jmax
    
      if(ind(j)==0) then
	ind(j) = 4
	soot(1, j) = j + 1 
	soot(2, j) = j + l 
	soot(3, j) = j - 1
	soot(4, j) = j - l
      endif
      
    enddo
    
    end subroutine get_soot_el
    subroutine get_JE(jmax, j, tau, je)
      implicit none
      integer, intent(in)		:: jmax
      complex*16, intent(in)	:: j(2*jmax)
      real*8, intent(in)		:: tau(3,2,jmax)
      complex*16, intent(out)	:: je(3,jmax)
      
      integer i, m
      m = 1
      do i = 1, jmax
	je(:,i)  = j(m) * tau(:,1,i) + j(m+1) * tau(:,2,i)
	m = m + 2
      enddo
    end subroutine get_JE 
    
end