ccccccc
        module potential
            implicit none
            contains                                
ccccccc
!V=V0/(1+exp(r-r0)/a)
            complex*16 function WSpot(v0,a,r0,r)
                real*8::v0,a,r0,r
                WSpot=-v0/(1d0+exp((r-r0)/a))
            end function WSpot
        end module