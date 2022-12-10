module gas

  use teio,       only: SPCSLEN,LINELEN

  use io,         only: io_adjustl,&
                        io_unit

  use input,      only: InputData, &
                        print_err
  
  use species,    only: find_spec_index_list,&
                        add_species

  use mesh,       only: get_neighb,&
                        Grid,&
                        Coor,&
                        SiteNum
                        
  use statistics, only: update_event_num

  implicit none
  private
  save

  public  :: init_gas,   &
             adsorption, &
             print_gas_info


  !private :: 

  type, public :: GasType
    character(len=SPCSLEN)    :: GName  !Gas name
    double precision          :: Prs   !Gas pressure
    double precision          :: AdsorbRate !Gas 
    double precision          :: Prob       ! 
    integer                   :: ProdNum
    integer,                dimension(2)  :: ProdI !adsorption species index
  end type GasType

  type(GasType), dimension(:), allocatable :: GasList
  integer                                  :: GasNum
  double precision, public                 :: GasProb

contains
  
  subroutine init_gas(inp) !initiate GasList
    type(InputData), intent(in)  :: inp
    integer          :: i

    GasProb=0
    GasNum=inp%GasNum
    if(GasNum==0) return

    allocate(GasList(inp%GasNum))
    do i=1,GasNum
      GasList(i)%GName=inp%GasName(i)
      GasList(i)%Prs=inp%GasPrs(i)
      GasList(i)%AdsorbRate=inp%GasAdsorbRate(i)
      call find_spec_index_list(inp%GasProd(i,2),GasList(i)%ProdI(1))
      if(trim(adjustl(inp%GasProd(i,3)))/='') then
        GasList(i)%ProdNum=2
        call find_spec_index_list(inp%GasProd(i,3),GasList(i)%ProdI(2))
      else
        GasList(i)%ProdNum=1
      end if
      GasList(i)%Prob=inp%GasPrs(i)*inp%GasAdsorbRate(i)*SiteNum
      GasProb=GasProb+GasList(i)%Prob
    end do

   end subroutine init_gas

   !-----------------------------------------------------------!

   subroutine adsorption()  ! for dense case, this subroutine needs revision.
     integer          :: id,i
     double precision :: probr,probs,randn
     call random_number(randn)
     probr=GasProb*randn
     probs=0
     do i=1,GasNum  
       probs=probs+GasList(i)%Prob
       if(probs>probr) then
         id=i
         goto 99
       end if
     end do
     call print_err("Error in choosing gas adsorption!")
99   do  i=1,GasList(id)%ProdNum
       !write(*,*) "adsorption gas"
       !write(*,*) GasList(id)%ProdI(i)
       call add_species(GasList(id)%ProdI(i))
     end do
     call update_event_num("adsorption",id)
   end subroutine adsorption

   !-----------------------------------------------------------!

   subroutine print_gas_info(u)
     implicit none
     integer,intent(in) :: u
     integer            :: i,j

     write(u,*) "GasNum:  ",io_adjustl(GasNum)
     write(u,*) "GasName/Pressure/AdsorbRate/AdsorbProduct"
     do i=1, GasNum
         if(GasList(i)%ProdNum==2) then
             write(u,'(" ",A5," ",E10.3," ",E11.4," ",I3," and",I3)') &
                GasList(i)%GName,GasList(i)%Prs,GasList(i)%AdsorbRate,&
                GasList(i)%ProdI(1),GasList(i)%ProdI(2)
         else if(GasList(i)%ProdNum==1) then
             write(u,'(" ",A5," ",E10.3," ",E11.4,"    ",I3)') &
                GasList(i)%GName,GasList(i)%Prs,GasList(i)%AdsorbRate,&
                GasList(i)%ProdI(1)
         end if
     end do
     write(u,*)

   end subroutine print_gas_info

end module gas
