ccccccc        
        program main
            use rmatmod 
            use input
            implicit none
ccccccc
            real*8::t1,t2
ccccccc
            call cpu_time(t1)
            call get_info()
ccccccc
            call readinput()  
ccccccc
            call getchannelalphaD()
ccccccc
            call rmat_int()
ccccccc
            call getpot(str)
ccccccc
            call rmatrix()
ccccccc
            deallocate(Vc)
            call cpu_time(t2)
ccccccc
            write(*,*) 'Running time: ',t2-t1
ccccccc
        end program 