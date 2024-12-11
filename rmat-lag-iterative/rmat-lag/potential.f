ccccccc
        module potential
            use parameter  
            implicit none
            contains                                
ccccccc
!V=V0/(1+exp(r-r0)/a)
            complex*16 function WSpot(v0,a,r0,r)
                real*8::v0,a,r0,r
                WSpot=-v0/(1d0+exp((r-r0)/a))
            end function WSpot
ccccccc
            real*8 function gausspot(r,v0,r0,a)
                real*8::r,v0,r0,a
                if (a.gt.1e-6) then
                gausspot=V0*exp(-(r-r0)**2/a**2)
                else
                write(*,*)'a too small in gausspot!'
                stop
                end if
            end function gausspot
ccccccc
            complex*16 function testpot(r)
                real*8::r
                testpot=-8d0*exp(-0.16d0*r**2)+4d0*exp(-0.04d0*r**2)
            end function
ccccccc

            function dws(r,v0,r0,a)
                implicit none
                real*8 r,v0,r0,a,dws,aux
                  if (r<1e-6) r=1e-6
                     if (a>1e-6) then
                       aux=exp(-(r-r0)/a)
                       dws=-v0*aux/(1d0+aux)**2/a
                     else
                       write(0,*)'derivative WS: a too small!';stop
                   endif
                   return
                end function

          
        end module
