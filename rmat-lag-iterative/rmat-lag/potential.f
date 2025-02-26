ccccccc
        module potential
            use parameter  
            implicit none
            contains                                
ccccccc
!V=V0/(1+exp(r-r0)/a)
            complex*16 function WSpot(v0,a0,r0,r)
                real*8::v0,a0,r0,r
                WSpot=-v0/(1d0+exp((r-r0)/a0))
            end function WSpot
ccccccc
            complex*16 function vpotn(v0,a0,r0,r)
                implicit none 
                real*8::v0,a0,r0,r
                vpotn=-v0*(1+cosh(r0/a0))/(cosh(r/a0)+cosh(r0/a0))
            end function vpotn
ccccccc
            complex*16 function vcoul(z12,r0,r)
                implicit none
                real*8::z12,r0,r
                if(r>=r0) then 
                vcoul=z12*e2/r
                else
                vcoul=z12*e2/2d0/r0*(3d0-r**2/r0**2)    
                end if
            end function vcoul
ccccccc
            complex*16 function vpot(v0,a0,r0,z12,r)
                implicit none 
                real*8::v0,a0,r0,z12,r
                vpot=vpotn(v0,a0,r0,r)+vcoul(z12,r0,r)
            end function vpot
        end module
