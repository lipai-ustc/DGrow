module core
   
  use input,   only: InputData, &
                     print_err, &
                     Kb

  use io,      only: io_adjustl,&
                     io_lower

  use mesh,    only: MeshSize,  &
                     CoreSize,  &
                     Grid,      &
                     period_y

  use statistics, only: update_event_num

  implicit none
  private 
  save

  public  :: init_core, &
             get_nn,    &
             att_core,  &
             det_core,  &
             save_core
                          
  private :: update_E,  &
             get_edge
                          
  type edge
      integer         :: enum    !edge sites
      integer         :: snum    !species num on edge
      integer         :: addnum
      integer         :: detnum
      integer         :: front   ! record the growth
      integer, dimension(:,:),allocatable   :: E ! i, (x,y,occ,nn,prob)
  end type

  type(edge), dimension(4), public :: Edges ! 1:4,left,right,bottom,top

  double precision, public :: EdgeProb
  double precision, dimension(:,:),allocatable :: DetProb
  double precision, dimension(:,:),allocatable :: MDetProb

  integer,dimension(2)          :: beg
  
contains

  !----------------------------------------------------------------!

  subroutine init_core(inp)
    implicit none
    type(InputData), intent(in)   :: inp
    integer                       :: i,j
    integer                       :: x,y

    allocate(DetProb(inp%SpecNum,3))
    allocate(MDetProb(inp%MSpecNum,3))
    DetProb=0
    MDetProb=0
    do i=1,inp%SpecNum
        do j=1,3
            if(inp%DetBrr(i,j)>0) then
                DetProb(i,j)=inp%Vib*exp(-inp%DetBrr(i,j)/Kb/inp%T)
            end if
        end do
    end do
    do i=1,inp%MSpecNum
        do j=1,3
            if(inp%MDetBrr(i,j)>0) then
                MDetProb(i,j)=inp%Vib*exp(-inp%MDetBrr(i,j)/Kb/inp%T)
            end if
        end do
    end do

    do i=1,4
        Edges(i)%snum=0
        Edges(i)%addnum=0
        Edges(i)%detnum=0
        Edges(i)%front=0
    end do

    if(period_y==0) then 
        if(inp%MeshType==1) then ! for hexagonal grid
            Edges(1)%enum=CoreSize(2)
            Edges(2)%enum=CoreSize(2)
            Edges(3)%enum=CoreSize(1)/2-1
            Edges(4)%enum=CoreSize(1)/2-1
            allocate(Edges(1)%E(CoreSize(2),5))
            allocate(Edges(2)%E(CoreSize(2),5))
            allocate(Edges(3)%E(CoreSize(1)/2-1,5))
            allocate(Edges(4)%E(CoreSize(1)/2-1,5))
            beg=MeshSize/2-CoreSize/2
            if(mod(beg(1),2)==0) beg(1)=beg(1)-1
            Edges(1)%E(:,1)=beg(1)
            Edges(2)%E(:,1)=beg(1)+CoreSize(1)-1
            do i=1,CoreSize(2)
                Edges(1)%E(i,2)=beg(2)+i-1
                Edges(2)%E(i,2)=beg(2)+i-1
            end do
            do i=1,CoreSize(1)/2-1
                Edges(3)%E(i,1)=beg(1)+2*i
                Edges(4)%E(i,1)=beg(1)+2*i-1
            end do
            Edges(3)%E(:,2)=beg(2)
            Edges(4)%E(:,2)=beg(2)+CoreSize(2)-1
        else if(inp%MeshType==2) then ! for square grid
            Edges(1)%enum=CoreSize(2)
            Edges(2)%enum=CoreSize(2)
            Edges(3)%enum=CoreSize(1)
            Edges(4)%enum=CoreSize(1)
            allocate(Edges(1)%E(CoreSize(2),5))
            allocate(Edges(2)%E(CoreSize(2),5))
            allocate(Edges(3)%E(CoreSize(1),5))
            allocate(Edges(4)%E(CoreSize(1),5))
            beg=MeshSize/2-CoreSize/2

            Edges(1)%E(:,1)=beg(1)
            Edges(2)%E(:,1)=beg(1)+CoreSize(1)-1
            do i=1,CoreSize(2)
                Edges(1)%E(i,2)=beg(2)+i-1
                Edges(2)%E(i,2)=beg(2)+i-1
            end do

            do i=1,CoreSize(1)
                Edges(3)%E(i,1)=beg(1)+i
                Edges(4)%E(i,1)=beg(1)+i-1
            end do
            Edges(3)%E(:,2)=beg(2)
            Edges(4)%E(:,2)=beg(2)+CoreSize(2)-1
        else 
            call print_err("invalid MeshType in init core")
        end if

        do i=1,4
            Edges(i)%E(:,3)=0 !occ , 0 for void, positive value for species, negative for mspecies
            Edges(i)%E(:,4)=1 !nn
            Edges(i)%E(:,5)=0 !prob
        end do

        do i=1,4
            do j=1,Edges(i)%enum
                Grid(Edges(i)%E(j,1),Edges(i)%E(j,2))%occ=-1 !-1 for available, -10 for attached 
            end do
        end do

        ! restore edge from Status file
        if(trim(io_lower(inp%Restart))=='true')  then 
            do i=1,4
                if(inp%EdgeSiteNum(i)/=Edges(i)%enum) then
                    call print_err("inp%EdgeSiteNum(i)/=Edges(i)%enum in init_core")
                end if
            end do
            do i=1,Edges(1)%enum 
                if(inp%Edge1(i)/=0) then
                    call att_core(inp%Edge1(i),Edges(1)%E(i,1),Edges(1)%E(i,2),0)
                end if
            end do
            do i=1,Edges(2)%enum 
                if(inp%Edge2(i)/=0) then
                    call att_core(inp%Edge2(i),Edges(2)%E(i,1),Edges(2)%E(i,2),0)
                end if
            end do
            do i=1,Edges(3)%enum 
                if(inp%Edge3(i)/=0) then
                    call att_core(inp%Edge3(i),Edges(3)%E(i,1),Edges(3)%E(i,2),0)
                end if
            end do
            do i=1,Edges(4)%enum 
                if(inp%Edge4(i)/=0) then
                    call att_core(inp%Edge4(i),Edges(4)%E(i,1),Edges(4)%E(i,2),0)
                end if
            end do
        end if
    else if(period_y==1) then 
        if(inp%MeshType==1.or.inp%MeshType==2) then
            Edges(1)%enum=CoreSize(2)
            allocate(Edges(1)%E(CoreSize(2),5))
            Edges(1)%E(:,1)=CoreSize(1)
            do i=1,MeshSize(2)
                Edges(1)%E(:,2)=i
            end do
        else 
            call print_err("error in init_core with invalid inp%MeshType")
        end if

        Edges(1)%E(:,3)=0 !occ , 0 for void, positive value for species, negative for mspecies
        Edges(1)%E(:,4)=1 ! which att
        Edges(1)%E(:,5)=0 !prob

        do i=1,Edges(1)%enum
            Grid(Edges(1)%E(i,1),Edges(i)%E(i,2))%occ=-1 !-1 for available, -10 for attached 
        end do
        
    else 
        call print_err("error in init_core with invalid period_y")
    end if
            
  end subroutine init_core

  !--------------------------------------------------------------------!

  subroutine refresh(edgei)
    implicit none
    integer,      intent(in)       :: edgei
    integer   :: i
    
    ! refresh Edges(edgei)%snum,front
    ! refresh Edges(edgei)%E( 3:occ, 4:nn, 5:prob)
    ! refresh Grid,EdgeProb
    Edges(edgei)%snum=0
    Edges(edgei)%front=Edges(edgei)%front+1

    do i=1,Edges(edgei)%enum
        Edges(edgei)%E(i,3)=0
        Grid(Edges(edgei)%E(i,1),Edges(edgei)%E(i,2))%occ=-1
        EdgeProb=EdgeProb-Edges(edgei)%E(i,5)
    end do

    do i=1,Edges(edgei)%enum
        call update_E(edgei,i) ! E( 4:nn, 5: prob) will be updated here
    end do
  
  end subroutine 
  
  !--------------------------------------------------------------------!

  subroutine get_nn(x,y,nn)
    implicit none
    integer,       intent(in)      :: x,y
    integer,       intent(out)      :: nn
    integer                        :: edgei,ind
    call get_edge(x,y,edgei,ind)
    nn=Edges(edgei)%E(ind,4)
  end subroutine get_nn

  !--------------------------------------------------------------------!

  subroutine update_E(edgei,ind)
    implicit none
    integer,       intent(in)      :: edgei,ind
    integer      :: nn,indnn,i
    nn=1
    if(period_y==0) then
        if(ind==0.or.ind==Edges(edgei)%enum+1) return
        if(ind==1) then
            nn=nn+1
        else if(Edges(edgei)%E(ind-1,3)/=0) then  !-10 for occupied
            nn=nn+1
        end if
        if(ind==Edges(edgei)%enum) then
            nn=nn+1
        else if(Edges(edgei)%E(ind+1,3)/=0) then
            nn=nn+1
        end if
    else if(period_y==1) then
        if(ind==0) indnn=Edges(edgei)%enum
        if(ind==Edges(edgei)%enum+1) indnn=1
        if(ind==1) then  !periodic condition
            indnn=Edges(edgei)%enum
        else
            indnn=ind-1
        end if
        if(Edges(edgei)%E(indnn,3)/=0) then  !-10 for occupied
            nn=nn+1
        end if

        if(ind==Edges(edgei)%enum) then
            indnn=1
        else 
            indnn=ind+1
        end if
        if(Edges(edgei)%E(indnn,3)/=0) then  !-10 for occupied
            nn=nn+1
        end if
    end if

    Edges(edgei)%E(ind,4)=nn

    if(Edges(edgei)%E(ind,3)==0) then
        Edges(edgei)%E(ind,5)=0
    else
        EdgeProb=EdgeProb-Edges(edgei)%E(ind,5)
        if(Edges(edgei)%E(ind,3)>0) then
            Edges(edgei)%E(ind,5)=DetProb(Edges(edgei)%E(ind,3),nn)
        else 
            Edges(edgei)%E(ind,5)=MDetProb(-Edges(edgei)%E(ind,3),nn)
        end if
        EdgeProb=EdgeProb+Edges(edgei)%E(ind,5)
    end if
        
  end subroutine update_E

  !--------------------------------------------------------------------!

  subroutine get_edge(x,y,edgei,ind)
    implicit none
    integer,       intent(in)      :: x,y
    integer,       intent(out)      :: edgei, ind
    integer  i,j
    if(period_y==0) then
        if(x==beg(1)) then   ! left
            edgei=1
            ind=y-beg(2)+1
        else if(x==beg(1)+CoreSize(1)-1) then  !right
            edgei=2
            ind=y-beg(2)+1
        else if(y==beg(2)) then  !bottom
            edgei=3
            ind=(x-beg(1))/2
        else if(y==beg(2)+CoreSize(2)-1) then !top
            edgei=4
            ind=(x-beg(1)+1)/2
        else
            write(*,*) x,y,Grid(x,y)%occ
            write(*,*) "beg",beg
            call print_err("Error in get_edge in core module")
        end if
        if(ind<1.or.ind>Edges(edgei)%enum) call print_err("Error in get_edge in core module")
    else if(period_y==1) then
        if(x/=CoreSize(1)) call print_err("Error in get_edge")
        edgei=1
        ind=y
    end if
  end subroutine get_edge

  !--------------------------------------------------------------------!

  subroutine att_core(id,x,y,flag)
    implicit none
    integer,       intent(in)      :: id
    integer,       intent(in)      :: x,y
    integer,optional, intent(in)   :: flag
    integer                        :: i,j,k
    integer                        :: edgei,ind
    call get_edge(x,y,edgei,ind)
    Edges(edgei)%E(ind,3)=id
    Grid(x,y)%occ=-10
    Edges(edgei)%snum=Edges(edgei)%snum+1
    Edges(edgei)%addnum=Edges(edgei)%addnum+1
    do i=-1,1
        call update_E(edgei,ind+i)
    end do
    if(present(flag)) return
    call update_event_num("attachment",id)

    ! if fully growthed 
    if(Edges(edgei)%snum==Edges(edgei)%enum) call refresh(edgei) 
  end subroutine att_core

  !----------------------------------------------------------------!

  subroutine det_core(edgei,ind,id)
    implicit none
    integer,       intent(in)      :: edgei,ind,id
    integer ::i
    Edges(edgei)%E(ind,3)=0
    Grid(Edges(edgei)%E(ind,1),Edges(edgei)%E(ind,2))%occ=-1
    EdgeProb=EdgeProb-Edges(edgei)%E(ind,5)
    do i=-1,1
        call update_E(edgei,ind+i)
    end do
    Edges(edgei)%snum=Edges(edgei)%snum-1
    Edges(edgei)%detnum=Edges(edgei)%detnum+1
    call update_event_num("detachment",id)
  end subroutine det_core

  !----------------------------------------------------------------!
  
  subroutine save_core(u)
    implicit none
    integer, intent(in)       :: u
    integer                   :: i,j
    character(5),dimension(4) :: chars=(/'Edge1','Edge2','Edge3','Edge4'/)

    write(u,'(A)') "EdgeSiteNum"
    do i=1,4
        write(u,'(4(I5," "))') Edges(:)%enum
    end do
    write(u,*)
    write(u,'(A)') "EdgeOccNum"
    do i=1,4
        write(u,'(4(I5," "))') Edges(:)%snum
    end do
    write(u,*)
    do i=1,4
        write(u,'(A)') chars(i)
        do j=1,Edges(i)%enum
            write(u,'(A," "$)') trim(io_adjustl(Edges(i)%E(j,3)))
        end do
        write(u,*)
    end do
  end subroutine save_core
 
end module core
