program simple_cg
!  use constants
  implicit none

  integer, parameter         :: np  = 100      ! length of input data
  integer, parameter         :: na  = 8019     ! length of sparse matrix
  integer, parameter         :: nd  = 400*400  ! length of image data

  real, dimension(nd)    :: db          ! data vector
  real, dimension(nd)    :: Nsig        ! diagonal of noise matrix N
  real, dimension(np)    :: bvec,xvec   ! b,x vector
  real, dimension(nd,np) :: Amat        ! A matrix
  real, dimension(np,np) :: Pmat        ! P matrix

  real ptmp
  integer i,j,x,y

  ! ... read data vector ...
  db = 0.
  open(44,file='eg5.data')
  read(44,*)
  x=0
  do i = 1,26666
     read(44,*) db(x+1:x+6)
     x = x+6
  enddo
  read(44,*) db(x+1:x+4)
  close(44)
!  call write_data(nd,db,'data_vector_eg5')

  ! ... read image matrix ...
  Amat = 0.
  open(44,file='eg5.matrix')
  read(44,*)
  do i = 1,8019 ! for eg3
     read(44,*) x,y,ptmp
     if ((x<0).or.(x>nd-1)) stop 'asdioghadga'
     if ((y<0).or.(y>np-1)) stop '99d8aghdg'
     Amat(x+1,y+1) = ptmp
!!     print *, x,y,ptmp
  enddo
  close(44)
!  open(44,file='data_matrix_eg5')
!  do i = 1,nd
!     write(44,'(999e17.8)') Pmat(i,:)
!  enddo
!  close(44)

  ! ... read noise covariance "matrix" ...
  Nsig = 0.
  open(44,file='eg5.sigma')
  read(44,*)
  x = 0
  do i= 1,26666
     read(44,*) Nsig(x+1:x+6)
     x = x+6
  enddo
  read(44,*) Nsig(x+1:x+4)
  close(44)

  ! ... create Pmat = A^T*Ninv*A ...
  do i = 1,np
     do j = 1,np
        Pmat(j,i) = sum(Amat(:,i)*Amat(:,j)/Nsig(:))
     enddo
  enddo

  ! ... create bvec = A^T*Ninv*db ...
  do i = 1,np
     bvec(i) = sum(Amat(:,i)*db(:)/Nsig(:))
  enddo

  ! ... solve Pmat*x=bvec with CG solver ...
  xvec = 1.
  call cg(np,Pmat,xvec,bvec)
  call write_data(np,xvec,'cg_out_eg5')

  ! ... compare with true result ...
  bvec = 0.
  open(44,file='eg5.src')
  read(44,*)
  x=0
  do i = 1,16
     read(44,*) bvec(x+1:x+6)
     x = x+6
  enddo
  read(44,*) bvec(x+1:x+4)
  close(44)
  bvec = abs(bvec-xvec)/abs(bvec)
  call write_data(np,bvec,'cg_relerr_eg5')
  

end program simple_cg


! ----------------------------------------------------------

subroutine write_data(n,x,fn)
!  use constants
  implicit none

  integer, intent(in)                 :: n
  real, dimension(n), intent(out)     :: x
  character(len=*), intent(in)        :: fn

  integer i

  open(47,file=TRIM(fn))
  do i = 1,n
     write(47,'(i9,255e17.8)') i, x(i)
  enddo
  close(47)

  return
end subroutine write_data

! ----------------------------------------------------------

subroutine cg(np,A,x,b)  ! Solve Ax=b with OrthoMin algo
!  use constants
  implicit none

  integer, intent(in)                :: np
  real, dimension(np,np), intent(in) :: A
  real, dimension(np), intent(in)    :: b
  real, dimension(np), intent(inout) :: x

  real eps

  real, dimension(np) :: xi,vi,ri
  real                :: di,dip,wi,psi,d0
  integer i

  ! ... set some accuracy ...
  eps = 1.d-6
  ! ... initial state ...
  i = 0
  xi  = x  ! initial guess
  ri  = b - MATMUL(A,xi)
  vi  = ri
  di  = DOT_PRODUCT(ri,ri)
  dip = DOT_PRODUCT(vi,MATMUL(A,vi))
!!  write(*,*) 'initial residual: ', i, di, dip

  ! ... CG MINRES iteration ...
  do
     wi  = di/dip
     ri  = ri - MATMUL(A,vi)*wi
     xi  = xi + vi*wi
     d0  = di
     di  = DOT_PRODUCT(ri,ri)
     psi = -di/d0
     vi  = ri-vi*psi
     dip = DOT_PRODUCT(vi,MATMUL(A,vi))
     i   = i+1
!!     write(*,*) 'current residual: ', i, di, dip
     if (di<eps) exit
  enddo

  ! ... return final result ...
  x = xi

  return
end subroutine cg
