ccccccc 
!this subroutine is used to calculate the inverse of a matrix
!with the subroutine of zgetrf and zgetri
!which comes from the lapack library
!for an arbitrary complex square matrix
ccccccc
            subroutine mat_inv(a,lda,n)
            integer::n,lda,info,lwork
            complex*16::a(n,n)
            complex*16,allocatable::work(:)  
            integer,allocatable::ipiv(:)

            lwork=n*n
            allocate(work(lwork),ipiv(n))

ccccccc
            call zgetrf(n, n, a, lda, ipiv, info)
ccccccc
            if (info .ne. 0) then
                print *, 'ZGETRF failed, info = ', info
                stop
            end if
ccccccc
            call zgetri(n, a, lda, ipiv, work, lwork, info)
ccccccc
            if (info .ne. 0) then
                print *, 'ZGETRI failed, info = ', info
                stop
            end if
ccccccc
            deallocate(work,ipiv)
ccccccc
            end subroutine mat_inv
ccccccc
!complex eigenvalue solutions from lapack        
!n:dimension of matrix,a:matrix
!w:eigenvalues,vl,vr: left & right eigenvectors
            subroutine ZGEEVS(n,a,w,vl,vr)             
                implicit none 
                integer::n
                integer::info,lda,ldvr,ldvl,lwork
                complex*16::a(n,n),w(n),vl(n,n),vr(n,n)
                complex*16::work(2*n),rwork(2*n)
            !    allocate(work(lwork),rwork(lwork))
                ldvr=n
                ldvl=n
                lda=n
                lwork=2*n
                call zgeev('N','V',n,A,lda,w,vl,ldvl,vr,ldvr,
     &       work,lwork,rwork,info)
               ! write(*,*) 'info=',info
               ! deallocate(work,rwork)
            end subroutine