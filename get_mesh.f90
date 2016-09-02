module get_mesh
implicit none
 contains

!! Процедура получения количества ячеек

  subroutine getjmax(filename, jmax)
  implicit none
    character (len=32), intent(in):: filename
    integer, intent(out):: jmax
    integer :: i, k, ios
    real*8 :: dummy
    open(15, FILE=filename) 
    do i = 1, 5
      read(15,*) dummy
    enddo
    k = 0
    ios = 0
    do while (ios /= -1)
      read(15,*,iostat = ios) dummy
      k = k + 1  
    enddo
    jmax = k - 1
    close(15)
  end subroutine getjmax
  
  subroutine getjmax_mid(filename, jmax, l, h)
  implicit none
    character (len=32), intent(in):: filename
    integer, intent(out):: jmax, l, h
    integer :: i
    real*8 :: dummy
    open(15, FILE=filename) 
    do i = 1, 3
      read(15,*) dummy
    enddo
    read(15,*) jmax
    read(15, *) h, l
    close(15)
  end subroutine getjmax_mid  
!! Процедура сбора координат вершин ячеек  

  subroutine getrut(filename, jmax, r)
  implicit none
    character (len=32), intent(in) :: filename
    integer, intent(in) :: jmax
    real*8, intent(out) :: r(3,4,jmax)
    integer :: i, k
    real*8 :: dummy
    open (15, FILE=filename)
    do i = 1, 5
      read(15,*) dummy
    enddo
    do i = 1, jmax
      read(15,*) r(:,:,i)
      r(:,:,i) = r(:,:,i)*10d0
    enddo
    close(15)
  end subroutine getrut
  
  subroutine get_rut12(filename, jmax, r1, r2)
  implicit none
    character (len=32), intent(in) :: filename
    integer, intent(in) :: jmax
    real*8, intent(out) :: r1(3,4,jmax),r2(3,4,jmax)
    integer :: i, k
    real*8 :: dummy
    open (15, FILE=filename)
    do i = 1, 5
      read(15,*) dummy
    enddo
    
    do i = 1, jmax
      read(15,*) r1(:,:,i)
      r1(:,:,i) = r1(:,:,i)!*10d0
    enddo
    
    do i = 1 , jmax
      read(15,*) r2(:,:,i)
      r2(:,:,i) = r2(:,:,i)!*10d0
    enddo
    close(15)
  end subroutine get_rut12
  
!! Процедура сбора координат точек коллокации  

  subroutine getrht(jmax,rut,rht)
  implicit none
    integer, intent(in) :: jmax
    real*8, intent(in) :: rut(3,4,jmax)
    real*8, intent(out) :: rht(3,jmax)
    integer :: i, j
    do j = 1, jmax
      rht(:,j) = 0d0
      do i = 1, 4
	rht(:,j) = rut(:,i,j) + rht(:,j)
      enddo
      rht(:,j) = rht(:,j)/4d0
    enddo
  end subroutine getrht
  
!! Процедура построения локального базиса  

  subroutine get_LOCAL_CF(jmax, rut, norm, tau, area)
  implicit none
    integer, intent(in) :: jmax
    real*8, intent(in) :: rut(3,4,jmax)
    real*8, intent(out) :: norm(3,jmax), tau(3,2,jmax), area(jmax)
    real*8 :: RS(3), PG(3), vectRSPG(3), dnrm2
    integer :: i, j
    
    do j = 1, jmax
    
      RS(:) =  ((rut(:,3,j) + rut(:,2,j)) - (rut(:,1,j) + rut(:,4,j)))/2d0
      PG(:) =  ((rut(:,3,j) + rut(:,4,j)) - (rut(:,1,j) + rut(:,2,j)))/2d0
      
      vectRSPG(1) = (RS(2) * PG(3) - RS(3) * PG(2))
      vectRSPG(2) = (RS(3) * PG(1) - RS(1) * PG(3)) 
      vectRSPG(3) = (RS(1) * PG(2) - RS(2) * PG(1))
      
      area(j) = sqrt(vectRSPG(1)**2 + vectRSPG(2)**2 + vectRSPG(3)**2)
      
      norm(:,j) = vectRSPG(:)/area(j)
      
      tau(:,1,j) = PG(:) / (sqrt(PG(1)**2 + PG(2)**2 + PG(3)**2))
      
      tau(1,2,j) = norm(2,j) * tau(3,1,j) - norm(3,j) * tau(2,1,j)
      tau(2,2,j) = norm(3,j) * tau(1,1,j) - norm(1,j) * tau(3,1,j)
      tau(3,2,j) = norm(1,j) * tau(2,1,j) - norm(2,j) * tau(1,1,j)
      
    enddo
  endsubroutine get_LOCAL_CF
  
!! Процедура поиска максимального диаметра разбиения  

  subroutine get_hmax(jmax, rut, hmax)
  implicit none
    integer, intent(in) :: jmax
    real*8, intent(in) :: rut(3,4,jmax)
    real*8, intent(out) :: hmax
    real*8:: mod_h, h(3)
    integer :: i, j
    hmax = 0d0
    
    do j = 1, jmax
      h = rut(:,4,j) - rut(:,1,j)
      mod_h = sqrt(h(1)*h(1) + h(2)*h(2) + h(3)*h(3))
      if (mod_h > hmax) then 
	hmax = mod_h
      endif
      do i = 1, 3
	h = rut(:,i+1,j) - rut(:,i,j)
	mod_h = sqrt(h(1)*h(1) + h(2)*h(2) + h(3)*h(3))
	if (mod_h > hmax) then 
	  hmax = mod_h
	endif
      enddo
    enddo 
  end subroutine get_hmax

  
  subroutine fix_ord_cells(rut, l, h, jmax)
    
    implicit none
    integer, intent(in) :: jmax, l, h
    real*8, intent(inout) :: rut(3,4,jmax)
    real*8 :: tmp(3), tmp2(3,4)
    integer :: i, m, k
    
    do k = 1, h
      do m = 1, l/2
	tmp2 = rut(:,:,(k-1)*l + m)
	rut(:,:,(k-1)*l + m) = rut(:,:,(k)*l - m + 1)
        rut(:,:,(k)*l - m + 1) = tmp2
      enddo
    enddo
    
  endsubroutine

      
  subroutine get_mid (rut1,rut2,jmax,mid)
    implicit none
    
    integer, intent(in)	:: jmax
    real*8, intent(in)	:: rut1(3,4,jmax), rut2(3,4,jmax)  
    real*8, intent(out)	:: mid(3,4,jmax)
    integer :: i, j
    do i = 1, jmax
      do j = 1, 4
	mid(:,j,i) = (rut1(:,j,i) + rut2(:,4-j+1,i))/2
      enddo
    enddo
    
    
  endsubroutine
    
  subroutine WTF_rut(matrix,jmax,filename)
    implicit none
    integer, intent(in) :: jmax
    character(len=9), intent(in) :: filename
    real*8, intent(in) :: matrix(3,4,jmax)
    integer i
    open (10,file=filename)
    write(10,*) 1
    write(10,*) 1
    write(10,*) 1
    write(10,*) jmax
    write(10,*) 1, jmax
    
    do i = 1, jmax
      write(10,*) matrix(:,:,i)
    enddo
    
    close(10)
    
  end subroutine WTF_rut
  
  subroutine get_norm(jmax, rut, norm)
  implicit none
    integer, intent(in) :: jmax
    real*8, intent(in) :: rut(3,4,jmax)
    real*8, intent(out) :: norm(3,jmax)
    real*8 :: RS(3), PG(3), vectRSPG(3),area
    integer :: j
    
    do j = 1, jmax
    
      RS(:) =  ((rut(:,3,j) + rut(:,2,j)) - (rut(:,1,j) + rut(:,4,j)))/2d0
      PG(:) =  ((rut(:,3,j) + rut(:,4,j)) - (rut(:,1,j) + rut(:,2,j)))/2d0
      
      vectRSPG(1) = (RS(2) * PG(3) - RS(3) * PG(2))
      vectRSPG(2) = (RS(3) * PG(1) - RS(1) * PG(3)) 
      vectRSPG(3) = (RS(1) * PG(2) - RS(2) * PG(1))
      
      area = sqrt(vectRSPG(1)**2 + vectRSPG(2)**2 + vectRSPG(3)**2)
      
      norm(:,j) = vectRSPG(:)/area     
    enddo
  endsubroutine get_norm 
end