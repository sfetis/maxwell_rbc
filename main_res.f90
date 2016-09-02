
program main

use operators_int
use get_mesh
use AV_sub
use get_aij

implicit none

integer 		:: jmax, i, ALLOC_ERR, j, n0, nrhs, INFO, l, h, p
character (len=32) 	:: filename,  outname
character (len=4)	:: postfix
real*8 			:: x(3),g(3), y(3), eps, hmax
real*8 			:: start, finish, rb, k, phi, phi1, dnrm2,nrm

complex*16, allocatable	:: Aij(:,:), RHS(:,:), JE(:,:), JM(:,:)

real*8, allocatable 	:: rut2(:,:,:),	rut0(:,:,:),	rut1(:,:,:)
real*8, allocatable 	:: norm2(:,:),	norm0(:,:),	norm1(:,:)
real*8, allocatable 	:: rkt0(:,:)
real*8, allocatable 	:: tau0(:,:,:)
real*8, allocatable 	:: area0(:)

real*8, allocatable 	:: RCS(:)
integer, allocatable 	:: IPIV(:), soot(:,:), ind(:)
real*8 			:: pi


print *, 'getting mesh parameters'
filename = '20_80_c_05.dat'
pi = 3.1415926d0
nrhs = 360
k = 0.628d0

call getjmax_mid(filename, jmax, l, h)
print *, 'jmax=', jmax

!!********************************************************************************************************

allocate(rut1(3,4,jmax),rut0(3,4,jmax),rut2(3,4,jmax), soot(4,jmax), ind(jmax))

call get_rut12(filename, jmax, rut1, rut2)

call fix_ord_cells(rut2, l, h, jmax)

call get_mid (rut1, rut2, jmax, rut0)

call get_soot(jmax, h, l, soot, ind)



allocate(norm2(3,jmax),	norm0(3,jmax),	norm1(3,jmax))
allocate(rkt0(3,jmax))
allocate(tau0(3,2,jmax))
allocate(area0(jmax))

allocate(Aij(4*jmax, 4*jmax), RHS(4*jmax, nrhs), RCS(nrhs), IPIV(4*jmax), JE(3,jmax), JM(3,jmax))

call getrht (jmax, rut0, rkt0)

call get_LOCAL_CF (jmax, rut0, norm0, tau0, area0)

call get_norm (jmax, rut1, norm1)
call get_norm (jmax, rut2, norm2)

call get_hmax (jmax, rut0, hmax)
print *, 'hmax=', hmax
print *, '* * * * * * *' 

Aij = 0d0
print *, 'Start matrix construction...'

    call get_Ematr	(rut0, rkt0, tau0, norm0, norm1, k, hmax, jmax, Aij(1:2*jmax,1:2*jmax), 10,&
			  E_1, kernel_00, kernel_2)    
			  
    call get_Bmatr	(rut0, rkt0, tau0, norm0, norm1, k, hmax, jmax, Aij(1:2*jmax,(2*jmax+1):4*jmax), 10,&
			  E_1, kernel_3)
			  
    call get_Ematr	(rut0, rkt0, tau0, norm0, norm2, k, hmax, jmax, Aij((2*jmax+1):4*jmax,1:2*jmax), 10,&
			  E_1, kernel_00, kernel_2)    
			  
    call get_Bmatr	(rut0, rkt0, tau0, norm0, norm2, k, hmax, jmax, Aij((2*jmax+1):4*jmax,(2*jmax+1):4*jmax), 10,&
			  E_1, kernel_3)
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call get_da_matr	(norm0, norm1, rut0, tau0, area0, jmax, Aij(1:2*jmax, 1:2*jmax), k, soot, ind, 0)
   
    call get_db_matr	(norm0, norm1, tau0, jmax, Aij(1:2*jmax,(2*jmax+1):4*jmax), 0)
     
    call get_da_matr	(norm0, norm2, rut0, tau0, area0, jmax, Aij(2*jmax+1:4*jmax, 1:2*jmax), k, soot, ind, 1) 
     
    call get_db_matr	(norm0, norm2, tau0, jmax, Aij(2*jmax+1:4*jmax, 2*jmax+1:4*jmax), 1)
     
print *, 'Done'    
print *, '* * * * * * *'     


n0 = 4 * jmax
rb = hmax/10d0
print *, 'Start matrix inversion...'
call ZGETRF(n0, n0, Aij, n0, IPIV, INFO)
print *, 'Done'
print *, '* * * * * * *'  
!!********************************************************************************************************

  open (50, file = '20_80_c_05_3ghz_H.txt')
  print *, 'Start calculate RCS...'
! !!********************************************************************************************************
! !! Решение системы и получение результатов для каждого угла облучения
  do i = 1, nrhs 

    phi = ((2 * pi)/nrhs) * (i - 1)

    phi1 = (phi - pi) * 180d0/pi
    
    call get_mid_RHSH(phi, jmax, rkt0, norm1, norm2, tau0, k, RHS(:,i))
    
    call ZGETRS('N', n0, 1, Aij, n0, IPIV, RHS(:,i), n0, INFO)
    
    call get_RCS(phi, jmax, RHS(:, i), k, tau0, rkt0, area0, RCS(i))
    
    write(50,*) phi1, Rcs(i)

  enddo
  close(50)
!!********************************************************************************************************
  
  print *, 'Done'
  print *, '* * * * * * *' 
  
deallocate(norm2,norm0,norm1,stat = ALLOC_ERR)

deallocate(rut1, rut2, rut0, stat= ALLOC_ERR)

deallocate(rkt0,stat = ALLOC_ERR)

deallocate(tau0,stat = ALLOC_ERR)

deallocate(area0,stat = ALLOC_ERR)
  
deallocate (Aij, soot, ind, stat = ALLOC_ERR)

deallocate ( RHS, RCS, IPIV, JE, JM,stat = ALLOC_ERR)

end