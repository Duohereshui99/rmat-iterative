ccccccc
!the module store the delta function for integer numbers
      module deltaf
        implicit none


        contains

        integer function delta(i,j)
        implicit none
        integer, intent(in) :: i,j
        if (i==j) then
          delta=1
        else
          delta=0
        end if
        
        end function delta


      end module deltaf