program cristallo
    
integer,parameter::kr=selected_real_kind(12)
integer :: NAtomilato,i,j,k,l,Nbase
real :: d
real, dimension(:,:), allocatable :: base
logical ::cubi,fcc

     
     write*, "Quanti vettori base?" !vettori di base
     read*, Nbase
     allocate(base(3,0:Nbase-1)
     
     
     write*, "Inserisci 3*",Nbase,"posizioni (x,y,z) per la base"
     do j=0, Nbase-1
        read*, base(1,j), base(2,j), base(3,j)
     enddo
     
  write(unit=*,fmt="(a)",advance="no")" Ti piacciono i cubi semplici?"
  read*, cubi
  if (cubi) then
     write(unit=*,fmt="(a)",advance="no") "Quanti siti per lato vuoi?"
     read*, NAtomilato
     write(unit=*, fmt="(a)", advance="no") "Che parametro reticolare?"
     read*, d
!     write(unit=20,fmt=*) NAtomilato
    do j=0, NAtomilato-1
       do k=0, NAtomilato-1
          do l=0, NAtomilato-1
             do j=0,Nbase-1
                write(unit=20,fmt=*) j*d+base(1,j), k*d+base(2,j), l*d+base(3,j)
             enddo
          end do
       end do
    end do
 else
     write(unit=*,fmt="(a)",advance="no") "Preferisci un facce centrate?"
     read*, fcc
     if (fcc) then
         write(unit=*, fmt="(a)", advance="no") "Quanto di parametro reticolare?"
         read*, d   
     write(unit=*,fmt="(a)",advance="no") "Quanti siti per lato vuoi (lunghezza cella ripetuta)?"
     read*, NAtomilato


     
     
    do j=0, NAtomilato-1
       do k=0, NAtomilato-1
          do l=0, NAtomilato-1
             do i =0,Nbase-1
                write(unit=20,fmt=*) j*d+base(1,j), k*d+base(2,j), l*d+base(3,j)
                write(unit=20,fmt=*) (j*d+d/2)+base(1,j), (k*d+d/2)+base(2,j), l*d+base(3,j)
                write(unit=20,fmt=*) j*d+base(1,j), (k*d+d/2)+base(2,j), (l*d+d/2)+base(3,j)
                write(unit=20,fmt=*) (j*d+d/2)+base(1,j), k*d+base(2,j), (l*d+d/2)+base(3,j)

             end do
          end do
          
       end do
    enddo
    


     
     else
     write(unit=*, fmt="(a)", advance="no") "Allora BCC. Quanto di parametro reticolare?"
     read*, d
     
     write(unit=*,fmt="(a)",advance="no") "Quanti atomi per lato vuoi?"
     read*, NAtomilato
    do j=0, NAtomilato-1
       do k=0, NAtomilato-1
          do l=0, NAtomilato-1
             do j=0, Nbase-1
                write(unit=20,fmt=*) j*d+base(1,j), k*d+base(2,j), l*d+base(3,j)
                write(unit=20,fmt=*) (j*d+d/2)+base(1,j), (k*d+d/2)+base(2,j), (l*d+d/2)+base(3,j)
             enddo
             
          end do
       end do
    end do
     
     
 end if
end if


 
  end program cristallo
  
    
       
          
