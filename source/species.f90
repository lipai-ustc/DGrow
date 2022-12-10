module species

  use teio,    only: SPCSLEN,LINELEN

  use input,   only: print_err, &
                     InputData, &
                     Kb

  use mesh,    only: get_neighb,&
                     get_nb_direct,&
                     Coor,      &
                     MeshSize,  &
                     SiteNum,   &
                     Grid,      &
                     flagfixden,&
                     FixDenArea

  use io,      only: io_adjustl,&
                     io_lower,  &
                     io_center

  use statistics, only: SpecI,  &
                        update_event_num

  use core,    only: att_core,  &
                     det_core,  &
                     get_nn,    &
                     Edges,     &
                     EdgeProb

  implicit none
  private
  save

  public  ::  init_species,        &
              fix_den,             &
              find_spec_index_list,&
              find_spec_index_all, &
              add_species,         &
              detach,              &
              species_step,        &
              mspecies_step,       &
              mmmerge_step,        &
              save_species,        &
              print_species_info,  &
              print_species_density,&
              print_species_name,&
              species_prob_update

  private ::  init_event,          &
              find_rand_void_site, &
              del_species,         &
              !init_q,              &
              !del_q,               &
              diffusion,           &
              decomp,              &
              combine,             &
              update_around,       &
              ispecies_step,       &
              ispecies_prob_update,  &
              iispecies_prob_update, &
              mspecies_prob_update, &
              mmmerge_prob_update

  !------------------------------species info-----------------------------!

  type, public  :: SpecType  ! for species without MFA
     !character(len=SPCSLEN)         :: SName
     integer                        :: ID
     integer                        :: Num
     type(SpecI), pointer           :: head
     double precision               :: ProbSum
     integer                        :: DenFix
     integer                        :: FixNum
  end type SpecType
  type(SpecType), dimension(:), allocatable           :: SpecList

  integer   :: fixtmp(1000,2)

  type, public  :: MSpecType  ! for species with MFA
     !character(len=SPCSLEN)         :: SName
     integer                        :: ID
     integer                        :: Num
     double precision, dimension(:),allocatable      :: ProbDecomp
     double precision               :: ProbDesorp
     integer                        :: DenFix
     double precision               :: ProbAtt
  end type MSpecType
  type(MSpecType), dimension(:), allocatable         :: MSpecList
  
  double precision, dimension(:),allocatable         :: MMergeProb

  integer                                            :: SpecNum
  character(len=SPCSLEN), dimension(:), allocatable  :: SpecName
  integer                                            :: MSpecNum
  character(len=SPCSLEN), dimension(:), allocatable  :: MSpecName
  integer                                            :: AllSpecNum

  double precision, public                           :: SpecProb
  double precision, public                           :: MSpecProb
  double precision, public                           :: MMMergeProb

  !------------------------------event info-----------------------------!

  type, public :: SpecE
    double precision                 :: DiffBrr
    double precision                 :: DesorpBrr
    double precision                 :: DiffProb
    double precision                 :: DesorpProbP
    double precision                 :: AttBrr(3)
    double precision                 :: DetBrr(3)
    double precision                 :: AttProb(3)
    double precision                 :: DetProb(3)
    ! For Att/Det, the first dimension is for different species
    ! and the second dimension is for different local environments on the edge
    integer                                          :: DecompNum
    integer,             dimension(:,:), allocatable :: DecompRI  !DecompReaction reaction Index
    double precision,    dimension(:),   allocatable :: DecompBrr
    double precision,    dimension(:),   allocatable :: DecompProbP
    
    integer                                          :: MergeNum
    integer,             dimension(:,:), allocatable :: MergeRI  !Merge with Mspecies only
    double precision,    dimension(:),   allocatable :: MergeBrr
    double precision,    dimension(:),   allocatable :: MergeProbP
  end type SpecE

  type, public :: MSpecE
    double precision              :: DesorpBrr
    double precision              :: DesorpProbP
    integer                                          :: DecompNum
    integer,             dimension(:,:), allocatable :: DecompRI  !DecompReaction reaction Index
    double precision,    dimension(:),   allocatable :: DecompBrr
    double precision,    dimension(:),   allocatable :: DecompProbP
    double precision              :: AttBrr(3)
    double precision              :: DetBrr(3)
    double precision              :: AttProb(3)
    double precision              :: DetProb(3)
  end type MSpecE

  type(SpecE),  dimension(:), allocatable,  public   :: SpecEvt
  type(MSpecE), dimension(:), allocatable,  public   :: MSpecEvt

  ! for species. both species and mspecies is included in the matrix
  integer,  dimension(:,:), allocatable              :: MergeM   !MergeMatrix
  integer,  dimension(:,:), allocatable              :: MergeMI   !reactionindex
  double precision, dimension(:,:), allocatable      :: MergeMBrr  !MergeBrrMatrix
  double precision, dimension(:,:), allocatable      :: MergeMProb  
  
  ! for mspecies. only merge between two mspecies are included in the array
  integer                                            :: MMergeNum
  integer,  dimension(:,:), allocatable              :: MMergeRI !Merge reactant and product index
  double precision, dimension(:), allocatable        :: MMergeBrr !MergeBrr
  double precision, dimension(:), allocatable        :: MMergeProbP  !MergeProb

  ! for Desorp
  integer, dimension(:), allocatable                 :: DesorpI

  integer  :: mattflag

  !------------------------Event info End-------------------------------!

contains

  !---------------------------------------------------------------------!

  subroutine init_species(inp)
    implicit none
    type(InputData), intent(in)    :: inp
    integer                        :: i,j,k
    integer                        :: id1,id2,id3

    SpecProb=0
    MSpecProb=0
    SpecNum=inp%SpecNum
    MSpecNum=inp%MSpecNum
    AllSpecNum=SpecNum+MSpecNum
    if(SpecNum+MSpecNum<=0) call print_err("No species in the model!")

    ! for Species
    if(SpecNum>0) then
        allocate(SpecName(SpecNum),SpecList(SpecNum))
        SpecName=inp%SpecName
        ! UpdataSpecList
        do i=1,SpecNum
          !SpecList(i)%SName=SpecName(i)
          SpecList(i)%ID=i
          SpecList(i)%Num=0
          SpecList(i)%head=>null()
          SpecList(i)%ProbSum=0
          SpecList(i)%DenFix=0 ! will be reassigned after initial the given density 
        end do
    endif

    ! for MSpecies
    if(MSpecNum>0)then
        allocate(MSpecName(MSpecNum),MSpecList(MSpecNum))
        MSpecName=inp%MSpecName
        !UpdateMSpecList
        do i=1,MSpecNum
          !MSpecList(i)%SName=MSpecName(i)
          MSpecList(i)%ID=-i
          MSpecList(i)%Num=0
          MSpecList(i)%DenFix=0 ! will be reassigned after initial the given density 
          MSpecList(i)%ProbAtt=0
        end do
    endif

    call init_event(inp)

    ! Restart, restore species density
    if(trim(io_lower(inp%Restart))=='true')  then ! restore mesh from Status file
        do i=1, inp%GridOccNum
                call add_species(inp%GridOcc(i,3),inp%GridOcc(i,1),inp%GridOcc(i,2))
        end do
        do i=1, inp%MSpecNumR
            do j=1, inp%NumMSpecies(i)
                call add_species(-i)
            end do
        end do
    else   ! new job, set initial dens
        do i=1,SpecNum
            if(inp%IniSpecDen(i)>0) then
                do j=1,int(inp%IniSpecDen(i)*SiteNum)
                  call add_species(i,update=.False.)
                end do
            else if(inp%IniSpecNum(i)>0) then
                do j=1,inp%IniSpecNum(i)
                  call add_species(i,update=.False.)
                end do
            end if
        end do
        do i=1,MSpecNum
            if(inp%MIniSpecDen(i)>0) then
                do j=1,int(inp%MIniSpecDen(i)*SiteNum)
                  call add_species(-i,update=.False.)
                end do
            else if(inp%MIniSpecNum(i)>0) then
                do j=1,inp%IniSpecNum(i)
                  call add_species(-i,update=.False.)
                end do
            end if
        end do
    endif

    ! set DenFix after possible initialization on species density
    do i=1,MSpecNum
        MSpecList(i)%DenFix=inp%MSpecDenFix(i)
    end do
    do i=1,SpecNum
        SpecList(i)%DenFix=inp%SpecDenFix(i)
        SpecList(i)%FixNum=int(real(FixDenArea(1)*FixDenArea(2)*SpecList(i)%Num)/SiteNum)
    end do


    call species_prob_update()
    call mspecies_prob_update()
    call mmmerge_prob_update()

    mattflag=inp%mattflag
  end subroutine init_species

  !--------------------fix spec num in specific area-------------------!
  subroutine fix_den()
    implicit none
    integer   :: i,j,x,y
    integer   :: n,ind
    double precision :: rand_tmp
    if(flagfixden==0) return
    do i=1,SpecNum
        if(SpecList(i)%DenFix==1) then 
            ! find the number of spec i in fixarea
            n=0
            do x=FixDenArea(1),MeshSize(1)
                do y=1,MeshSize(2)
                    if(Grid(x,y)%occ==i) then
                        n=n+1
                        if(n<1000) then 
                            fixtmp(n,1)=x
                            fixtmp(n,2)=y
                        end if
                    end if
                end do
            end do
            if(n<SpecList(i)%FixNum) then
                do j=1,SpecList(i)%FixNum-n
                    call add_species(i)
                end do
            else if(n>SpecList(i)%FixNum) then
                do j=1,n-SpecList(i)%FixNum
                    call random_number(rand_tmp)
                    ind=ceiling(n*rand_tmp)
                    if(ind>1000) ind=ceiling(1000*rand_tmp)
                    call del_species(i,fixtmp(ind,1),fixtmp(ind,2))
                end do
            end if
        end if
    end do
  end subroutine fix_den

  !----------------------------Event info part--------------------------!
  subroutine  init_event(inp)
    implicit none
    type(InputData), intent(in)    :: inp

    integer                        :: i,j,k,id
    integer                        :: id1,id2,id3

    ! Updata Merge Matrix
    allocate(MergeM(AllSpecNum,AllSpecNum))
    allocate(MergeMI(AllSpecNum,AllSpecNum))
    allocate(MergeMBrr(AllSpecNum,AllSpecNum))
    allocate(MergeMProb(AllSpecNum,AllSpecNum))
    MergeM=0 ! 0 means no this merge reaction or the associative desorption event
    MergeMI=0 ! 0 means no this merge reaction or the associative desorption event
    MergeMBrr=-1.0  ! -1 means no this merge reaction
    MergeMProb=-1.0 ! -1 means no this merge reaction
    MMergeNum=0  ! count MMergeNum from 0
    do i=1,inp%MergeNum
      call find_spec_index_all(inp%MergeR(i,1),id1)
      call find_spec_index_all(inp%MergeR(i,2),id2)
      call find_spec_index_list(inp%MergeR(i,3),id3)
      MergeM(id1,id2)=id3
      MergeM(id2,id1)=id3
      MergeMI(id1,id2)=i
      MergeMI(id2,id1)=i
      MergeMBrr(id1,id2)=inp%MergeBrr(i)
      MergeMBrr(id2,id1)=inp%MergeBrr(i)
      if(id1<inp%SpecNum.and.id2<inp%SpecNum) then
         MergeMProb(id1,id2)=inp%Vib*exp(-inp%MergeBrr(i)/Kb/inp%T)*0.5
      else if(id1>inp%SpecNum.and.id2>inp%SpecNum) then
         MMergeNum=MMergeNum+1
         MergeMProb(id1,id2)=inp%Vib*exp(-inp%MergeBrr(i)/Kb/inp%T)
      else
         MergeMProb(id1,id2)=inp%Vib*exp(-inp%MergeBrr(i)/Kb/inp%T)
      end if
      MergeMProb(id2,id1)=MergeMProb(id1,id2)
    end do

    ! Update MMerge arrays
    if(MMergeNum>0) then
        allocate(MMergeRI(MMergeNum,4))
        allocate(MMergeBrr(MMergeNum))
        allocate(MMergeProbP(MMergeNum))
        allocate(MMergeProb(MMergeNum))
        k=1
        do i=1,inp%MergeNum
            call find_spec_index_list(inp%MergeR(i,1),id1)
            call find_spec_index_list(inp%MergeR(i,2),id2)
            call find_spec_index_list(inp%MergeR(i,3),id3)
            if(id1<0.and.id2<0) then
                MMergeRI(k,1)=id1
                MMergeRI(k,2)=id2
                MMergeRI(k,3)=id3
                MMergeRI(k,4)=i
                MMergeBrr(k)=inp%MergeBrr(i)
                MMergeProbP(k)=inp%Vib*exp(-inp%MergeBrr(i)/Kb/inp%T)
                k=k+1
            end if
        end do
    end if

    if(SpecNum>0) then
        allocate(SpecEvt(SpecNum))
        !UpdataSpecEvt
        do i=1,SpecNum
          ! Diffusion
          SpecEvt(i)%DiffBrr=inp%DiffBrr(i)
          if(inp%DiffBrr(i)>0) then
            SpecEvt(i)%DiffProb=inp%Vib*exp(-inp%DiffBrr(i)/Kb/inp%T)
          else
            SpecEvt(i)%DiffProb=-1
          end if
          SpecEvt(i)%DesorpBrr=-1   ! Desorption info will be updated later,
          SpecEvt(i)%DesorpProbP=-1  ! Assignment here for species without desorption event

          ! Attachment & Detachment
          if(inp%FlagCore) then
            SpecEvt(i)%AttBrr(:)=inp%AttBrr(i,:)
            SpecEvt(i)%DetBrr(:)=inp%DetBrr(i,:)
          else 
            SpecEvt(i)%AttBrr(:)=-1  ! negative barrier indicates impossible events
            SpecEvt(i)%DetBrr(:)=-1
          end if
          do j=1,3
              if(SpecEvt(i)%AttBrr(j)>0) then
                SpecEvt(i)%AttProb(j)=inp%Vib*exp(-inp%AttBrr(i,j)/Kb/inp%T)
              else 
                SpecEvt(i)%AttProb(j)=-1
              end if
              if(SpecEvt(i)%DetBrr(j)>0) then
                SpecEvt(i)%DetProb(j)=inp%Vib*exp(-inp%DetBrr(i,j)/Kb/inp%T)
              else 
                SpecEvt(i)%DetProb(j)=-1
              end if
          end do 

          ! Decomp
          SpecEvt(i)%DecompNum=0  ! count DecompNum
          do j=1,inp%DecompNum
            if(SpecName(i)==inp%DecompR(j,1))  SpecEvt(i)%DecompNum=SpecEvt(i)%DecompNum+1
          end do
          if(SpecEvt(i)%DecompNum>0) then
              allocate(SpecEvt(i)%DecompRI(SpecEvt(i)%DecompNum,3))
              allocate(SpecEvt(i)%DecompBrr(SpecEvt(i)%DecompNum))
              allocate(SpecEvt(i)%DecompProbP(SpecEvt(i)%DecompNum))
              k=1
              do j=1,inp%DecompNum
                if(SpecName(i)==inp%DecompR(j,1)) then
                  call find_spec_index_list(inp%DecompR(j,2),id1)
                  SpecEvt(i)%DecompRI(k,1)=id1
                  call find_spec_index_list(inp%DecompR(j,3),id2)
                  SpecEvt(i)%DecompRI(k,2)=id2
                  SpecEvt(i)%DecompBrr(k)=inp%DecompBrr(j)
                  SpecEvt(i)%DecompProbP(k)=inp%Vib*exp(-inp%DecompBrr(j)/Kb/inp%T)
                  SpecEvt(i)%DecompRI(k,3)=j
                  k=k+1
                end if
              end do
          end if

          ! Merge with MSpecies only (between this species and any MSpecies)
          SpecEvt(i)%MergeNum=0   !count MergeNum
          do j=1,inp%MergeNum
              call find_spec_index_list(inp%MergeR(j,1),id1)
              call find_spec_index_list(inp%MergeR(j,2),id2)
              if(id1==i) then
                  if(id2<0) then
                      SpecEvt(i)%MergeNum=SpecEvt(i)%MergeNum+1
                  end if
              end if
              if(id2==i) then
                  if(id1<0) then
                      SpecEvt(i)%MergeNum=SpecEvt(i)%MergeNum+1
                  end if
              end if
          end do
          allocate(SpecEvt(i)%MergeRI(SpecEvt(i)%MergeNum,3))
          allocate(SpecEvt(i)%MergeBrr(SpecEvt(i)%MergeNum))
          allocate(SpecEvt(i)%MergeProbP(SpecEvt(i)%MergeNum))
          k=1
          if(SpecEvt(i)%MergeNum>0) then
              do j=1,inp%MergeNum
                  call find_spec_index_list(inp%MergeR(j,1),id1)
                  call find_spec_index_list(inp%MergeR(j,2),id2)
                  call find_spec_index_list(inp%MergeR(j,3),id3)
                  if(id1==i) then
                      if(id2<0) then
                          SpecEvt(i)%MergeRI(k,1)=id2
                          SpecEvt(i)%MergeRI(k,2)=id3
                          SpecEvt(i)%MergeRI(k,3)=j
                          SpecEvt(i)%MergeBrr(k)=inp%MergeBrr(j)
                          if(inp%MergeBrr(j)<0) call print_err("inp%MergeBrr(j)<0")
                          SpecEvt(i)%MergeProbP(k)=inp%Vib*exp(-inp%MergeBrr(j)/Kb/inp%T)
                          k=k+1
                      end if
                  end if
                  if(id2==i) then
                      if(id1<0) then
                          SpecEvt(i)%MergeRI(k,1)=id1
                          SpecEvt(i)%MergeRI(k,2)=id3
                          SpecEvt(i)%MergeRI(k,3)=j
                          SpecEvt(i)%MergeBrr(k)=inp%MergeBrr(j)
                          if(inp%MergeBrr(j)<0) call print_err("inp%MergeBrr(j)<0")
                          SpecEvt(i)%MergeProbP(k)=inp%Vib*exp(-inp%MergeBrr(j)/Kb/inp%T)
                          k=k+1
                      end if
                  end if
              end do
          end if
        end do
    end if   !Spec info update end

    if(inp%MSpecNum>0) then
        allocate(MSpecEvt(inp%MSpecNum))
        ! UpdateMSpecEvt
        do i=1, inp%MSpecNum
            MSpecEvt(i)%DesorpBrr=-1   ! Desorption info will be updated later,
            MSpecEvt(i)%DesorpProbP=-1  ! Assignment here for species without desorption event
            ! Decomp
            MSpecEvt(i)%DecompNum=0
            do j=1,inp%DecompNum
                if(MSpecName(i)==inp%DecompR(j,1)) MSpecEvt(i)%DecompNum=MSpecEvt(i)%DecompNum+1
            end do
            if(inp%DecompNum>0) then
                allocate(MSpecEvt(i)%DecompRI(MSpecEvt(i)%DecompNum,3))
                allocate(MSpecEvt(i)%DecompBrr(MSpecEvt(i)%DecompNum))
                allocate(MSpecEvt(i)%DecompProbP(SpecEvt(i)%DecompNum))
                k=1
                do j=1,inp%DecompNum
                    if(MSpecName(i)==inp%DecompR(j,1)) then
                        call find_spec_index_list(inp%DecompR(j,2),id1)
                        MSpecEvt(i)%DecompRI(k,1)=id1
                        call find_spec_index_list(inp%DecompR(j,3),id2)
                        MSpecEvt(i)%DecompRI(k,2)=id2
                        MSpecEvt(i)%DecompBrr(k)=inp%DecompBrr(j)
                        if(inp%DecompBrr(j)<0) call print_err("inp%DecompBrr(j)<0")
                        MSpecEvt(i)%DecompProbP(k)=inp%Vib*exp(-inp%DecompBrr(j)/Kb/inp%T)
                        MSpecEvt(i)%DecompRI(k,3)=j
                        k=k+1
                    end if
                end do
            end if
            ! core part
            if(inp%FlagCore) then
                MSpecEvt(i)%AttBrr(:)=inp%MAttBrr(i,:)
                MSpecEvt(i)%DetBrr(:)=inp%MDetBrr(i,:)
            else 
                MSpecEvt(i)%AttBrr(:)=-1  ! negative barrier indicates impossible events
                MSpecEvt(i)%DetBrr(:)=-1
            end if
            do j=1,3
                if(MSpecEvt(i)%AttBrr(j)>0) then
                    MSpecEvt(i)%AttProb(j)=inp%Vib*exp(-inp%MAttBrr(i,j)/Kb/inp%T)
                    write(*,*) "MAttBrr, i,j, attprob",i,j,MSpecEvt(i)%AttProb(j)
                else 
                    MSpecEvt(i)%AttProb(j)=0
                end if
                if(MSpecEvt(i)%DetBrr(j)>0) then
                    MSpecEvt(i)%DetProb(j)=inp%Vib*exp(-inp%MDetBrr(i,j)/Kb/inp%T)
                    write(*,*) "MDetBrr, i,j, detprob",i,j,MSpecEvt(i)%DetProb(j)
                else 
                    MSpecEvt(i)%DetProb(j)=0
                end if
            end do 
        end do
    end if

    ! Desorp events for both Spec and MSpec
    allocate(DesorpI(AllSpecNum))
    DesorpI=0
    do i=1,inp%DesorpNum
      if(inp%DesorpBrr(i)<0) call print_err("inp%DesorpBrr(i)<0")
      call find_spec_index_list(inp%DesorpR(i),id)
      if(id>0) then
        SpecEvt(id)%DesorpBrr=inp%DesorpBrr(i)
        SpecEvt(id)%DesorpProbP=inp%Vib*exp(-inp%DesorpBrr(i)/Kb/inp%T)
        DesorpI(id)=i
      else if(id<0) then
        MSpecEvt(-id)%DesorpBrr=inp%DesorpBrr(i)
        MSpecEvt(-id)%DesorpProbP=inp%Vib*exp(-inp%DesorpBrr(i)/Kb/inp%T)
        DesorpI(SpecNum-id)=i
      end if
    end do

  end subroutine init_event

  !----------------------------Event info end---------------------------!

  subroutine find_spec_index_list(namei, id)
    character(len=SPCSLEN), intent(in)    :: namei
    integer,                intent(out)   :: id 
    integer                               :: i
    do i=1,SpecNum
      if(trim(adjustl(namei))==trim(adjustl(SpecName(i)))) then
          id=i
          return
      end if
    enddo
    do i=1,MSpecNum
      if(trim(adjustl(namei))==trim(adjustl(MSpecName(i)))) then
          id=-i
          return
      end if
    enddo
    id=0
    ! call print_err("Error in find_spec_index!",1)
    ! call print_err("Out of the range of possible species:"//namei,3)
  end subroutine find_spec_index_list

  !---------------------------------------------------------------------!

  subroutine find_spec_index_all(namei, id)
    character(len=SPCSLEN), intent(in)    :: namei
    integer,                intent(out)   :: id 
    call find_spec_index_list(namei,id)
    if(id<0) id=-id+SpecNum
  end subroutine find_spec_index_all

  !---------------------------------------------------------------------!

  subroutine add_species(id,x,y,update,ffix2)
    implicit none
    integer,            intent(in) :: id   !serial number of species type
    integer,  optional, intent(in) :: x   !coordination of this species
    integer,  optional, intent(in) :: y   !coordination of this species
    logical,  optional, intent(in) :: update  
    logical,  optional, intent(in) :: ffix2
    integer                        :: x_temp,y_temp  

    type(SpecI),    pointer        :: q

    if(id==0) return   ! corresponds to associative desorption
     
    if(id<0)  then  ! species with MFA
      if(MSpecList(-id)%DenFix==0) then ! 
          MSpecList(-id)%Num=MSpecList(-id)%Num+1
          call mspecies_prob_update(id)
          call mmmerge_prob_update()
      end if
      return
    end if

    ! species without MFA
    ! if coordinate is not given, species will be add at randion void site
    if(.not. present(x))  then
        call find_rand_void_site(x_temp,y_temp)
    else
        x_temp=x
        y_temp=y
    end if

    !init q
    !call init_q(id,x_temp,y_temp,q)
    allocate(q)
    q%SID=id
    q%Coor(1)=x_temp
    q%Coor(2)=y_temp
    allocate(q%ProbNB(Coor%NBNum))
    if(SpecEvt(id)%MergeNum>0) allocate(q%ProbMerge(SpecEvt(id)%MergeNum))
    if(SpecEvt(id)%DecompNum>0) allocate(q%ProbDecomp(SpecEvt(id)%DecompNum))
    q%ProbDesorp=SpecEvt(id)%DesorpProbP ! it might be negative
    q%ProbSum=0

    ! add q at head of the list
    q%prev=>null()
    if(associated(SpecList(id)%head))  then ! already has this type of species
        q%next=>SpecList(id)%head
        SpecList(id)%head%prev=>q
    else       ! new one
        q%next=>null()
    endif
    SpecList(id)%head=>q
    SpecList(id)%Num=SpecList(id)%Num+1

    Grid(x_temp,y_temp)%occ=id
    Grid(x_temp,y_temp)%spec=>q
    q=>null()

    if(present(update).and.(update.eqv..false.)) return

    ! update events around this q (including this q)
    call ispecies_prob_update(Grid(x_temp,y_temp)%spec)
    call update_around(x_temp,y_temp)
    
    ! if spec den fixed, del one spec after the add_species
    if(present(ffix2)) return
    if(SpecList(id)%DenFix/=0) call del_species(id=id,ffix1=.False.)

  end subroutine add_species

  !---------------------------------------------------------------------!

  !subroutine init_q(id_q,x_q,y_q,q)
  !  implicit none
  !  integer,  intent(in) :: id_q
  !  integer,  intent(in) :: x_q,y_q
  !  type(SpecI), pointer, intent(out) :: q
  !  allocate(q)
  !  q%SID=id_q
  !  q%Coor(1)=x_q
  !  q%Coor(2)=y_q
  !  allocate(q%ProbNB(Coor%NBNum))
  !  if(SpecEvt(id_q)%MergeNum>0) allocate(q%ProbMerge(SpecEvt(id_q)%MergeNum))
  !  if(SpecEvt(id_q)%DecompNum>0) allocate(q%ProbDecomp(SpecEvt(id_q)%DecompNum))
  !  q%ProbDesorp=SpecEvt(id_q)%DesorpProbP ! it might be negative
  !  q%ProbSum=0
  !end subroutine init_q


  !---------------------------------------------------------------------!

  subroutine find_rand_void_site(x,y)
    implicit none
    integer, intent(out)           :: x,y
    double precision               :: randn
    integer                        :: i

    ! This site is void, and it has at least one void neighbour site
    do  !while(.true.)
      call random_number(randn)
      x=ceiling(MeshSize(1)*randn)
      call random_number(randn)
      y=ceiling(MeshSize(2)*randn)
      if(Grid(x,y)%occ/=0) cycle
      call get_neighb(x,y)
      do i=1,Coor%NBNum
        if(Grid(Coor%NB(i,1),Coor%NB(i,2))%occ==0) return
      end do
    end do
  end subroutine find_rand_void_site

  !---------------------------------------------------------------------!

  subroutine del_species(id,x1,y1,update,ffix1)
    implicit none
    integer,            intent(in) :: id   !serial number of species type, index
    integer,  optional, intent(in) :: x1,y1   !coordination of this species
    logical,  optional, intent(in) :: update  ! update or not
    logical,  optional, intent(in) :: ffix1    ! fix flag
    integer                        :: x,y,n,i
    double precision               :: randn
    type(SpecI),    pointer                     :: p,q

    if(id<0)  then  ! species with MFA
      if(MSpecList(-id)%DenFix==0) then ! 
          MSpecList(-id)%Num=MSpecList(-id)%Num-1
          call mspecies_prob_update(id)
          call mmmerge_prob_update()
      end if
      return
    end if

    if(.not. present(x1))  then
       call random_number(randn)
       n=ceiling(SpecList(id)%Num*randn)
       p=>Speclist(id)%head
       do i=1,n-1
         p=>p%next
       end do
       x=p%Coor(1)
       y=p%Coor(2)
       p=>null()
    else
      x=x1
      y=y1
    end if

    if(id/=Grid(x,y)%occ) then
      call print_err("Error in del_species!",1)
      call print_err("id is: "//io_adjustl(id),2)
      call print_err("but Grid(x,y)%occ is: "//io_adjustl(Grid(x,y)%occ),2)
    end if


    if(Grid(x,y)%occ==0) then  ! revision is needed here
      call print_err("Error in del_spec in species module",1)
      call print_err("grid(x,y)%occ/=1",2)
      call print_err("No species at "//io_adjustl(x)//","//io_adjustl(y)//" when delete this species",3)
    endif

    ! delete species q
    q=>Grid(x,y)%spec
    SpecList(q%SID)%ProbSum=SpecList(q%SID)%ProbSum-q%ProbSum
    SpecProb=SpecProb-q%ProbSum
    Grid(x,y)%spec=>null()
    Grid(x,y)%occ=0
    SpecList(id)%Num=SpecList(id)%Num-1

    if(associated(q%prev)) then   ! q is not the first one in the chain
        if(associated(q%next)) then    ! q is not the last one in the chain 
            q%prev%next=>q%next
            q%next%prev=>q%prev
        else                           ! q is the last one in the chain
            q%prev%next=>null()
        end if
    else                          ! q is the first one
        if(associated(q%next)) then    ! q is not the last one in the chain 
            SpecList(q%SID)%head=>q%next
            q%next%prev=>null()
        else                           ! q is the last one in the chain
            SpecList(q%SID)%head=>null()
        end if
    end if


    ! dell q
    deallocate(q%ProbNB)
    if(allocated(q%ProbMerge)) deallocate(q%ProbMerge)
    if(allocated(q%ProbDecomp)) deallocate(q%ProbDecomp)
    !write(*,*) "in del_spec x,y,id: " , x,y,id
    deallocate(q)

    !write(*,*) "out del_spec x,y,id: " , x,y,id

    if(present(update).and.(update.eqv..false.))  return

    ! update events around this q
    call update_around(x,y)

    ! if spec den fixed, del one spec after the add_species
    if(present(ffix1)) return
    if(SpecList(id)%DenFix/=0) call add_species(id=id,ffix2=.False.)
    
  end subroutine del_species

  !---------------------------------------------------------------------!

  !subroutine del_q(q)
  !  implicit none
  !  type(SpecI), pointer, intent(out) :: q
  !  type(SpecI), pointer              :: p
  !  p=>q
  !  q=>null()
  !  deallocate(p%ProbNB)
  !  if(allocated(p%ProbMerge)) deallocate(p%ProbMerge)
  !  if(allocated(p%ProbDecomp)) deallocate(p%ProbDecomp)
  !  !deallocate(p)
  !end subroutine 

  !---------------------------------------------------------------------!

  subroutine diffusion(id,x,y,drct)
    implicit none
    integer,            intent(in) :: id    !serial number of species type, index
    integer,            intent(in) :: x,y   !coordination of this species
    integer,            intent(in) :: drct  !coordination of this species
    integer                        :: x1,y1

    call get_nb_direct(x,y,drct)
    x1=Coor%Drct(1)
    y1=Coor%Drct(2)

    !if(id/=Grid(x,y)%occ) then
    !  call print_err("Error in diffusion!",1)
    !  call print_err("id is: "//io_adjustl(id),2)
    !  call print_err("but Grid(x,y)%occ is: "//io_adjustl(Grid(x,y)%occ),2)
    !end if

    !if(Grid(x1,y1)%occ/=0) then
    !  call print_err("Error in diffusion!",1)
    !  call print_err("The accept site is not void. ",3)
    !end if 

    ! do not need to alloc or dealloc q for diffusion
    Grid(x1,y1)%occ=Grid(x,y)%occ
    Grid(x1,y1)%spec=>Grid(x,y)%spec
    Grid(x1,y1)%spec%Coor=Coor%Drct
    Grid(x,y)%occ=0
    Grid(x,y)%spec=>null()

    ! update events around this q
    call update_around(x,y)
    call update_around(x1,y1)

    call update_event_num("diffusion",id)

  end subroutine diffusion

  !---------------------------------------------------------------------!

  subroutine decomp(id,evt,x,y)
    implicit none
    integer,            intent(in) :: id    !serial number of species type, index
    integer,            intent(in) :: evt   !serial number of species type, index
    integer, optional,  intent(in) :: x,y   !coordination of this species
    integer                        :: drct,i

    if(id<0) then   ! MSpecies
      call del_species(id)
      call add_species(MSpecEvt(-id)%DecompRI(evt,1))
      call add_species(MSpecEvt(-id)%DecompRI(evt,2))
      call update_event_num("decomposition",MSpecEvt(-id)%DecompRI(evt,3))
      return
    end if
       
    call get_neighb(x,y)
    do i=1,Coor%NBNum
      if(Grid(Coor%NB(i,1),Coor%NB(i,2))%occ==0) then
        drct=i
        exit
      end if
      if(i==Coor%NBNum) then
        call print_err("Error in decomp!no place for the new species to occupy.")
      end if
    end do

    call del_species(id,x,y,update=.false.)
    call add_species(SpecEvt(id)%DecompRI(evt,1),x,y)
    call add_species(SpecEvt(id)%DecompRI(evt,2),Coor%NB(drct,1),Coor%NB(drct,2))

    call update_event_num("decomposition",SpecEvt(id)%DecompRI(evt,3))

  end subroutine decomp

  !---------------------------------------------------------------------!

  subroutine combine(id,x,y,drct,evt) ! merge is a keyword in fortran 
    implicit none
    integer,            intent(in) :: id        ! serial number of species type, index
    integer,            intent(in) :: x,y       ! coordination of this species
    integer, optional,  intent(in) :: drct      ! direct as index of another reactant
    integer, optional,  intent(in) :: evt       ! merge with mspecies
    integer                        :: id1,id2   !id of reactant2 and product1
    integer                        :: evt_id

    call del_species(id,x,y)
    if(present(drct)) then  ! merge with species
        ! get_nb_direct is needed here because Drct changed in update in del_species 
        call get_nb_direct(x,y,drct)
        if(Grid(Coor%Drct(1),Coor%Drct(2))%occ==0) &
            call print_err("Error in merge! Some reactant is void.")
        id1=Grid(Coor%Drct(1),Coor%Drct(2))%occ
        id2=MergeM(id,id1)
        evt_id=MergeMI(id,id1)
        call del_species(id1,Coor%Drct(1),Coor%Drct(2))
    else if(present(evt)) then  ! merge with mspecies
        id1=SpecEvt(id)%MergeRI(evt,1)
        id2=SpecEvt(id)%MergeRI(evt,2)
        evt_id=SpecEvt(id)%MergeRI(evt,3)
        call del_species(id1)
    end if
    if(id2<0) then
        call add_species(id2)     ! MSpecies
    else if(id2>0) then
        call add_species(id2,x,y) ! Species
    end if
    call update_event_num("merge",evt_id)

  end subroutine combine

  !---------------------------------------------------------------------!
  subroutine detach()
    implicit none
    integer :: i,j,k
    integer :: x,y,id
    double precision :: eprobr, eprobs, erandn
    call random_number(erandn)
    eprobr=EdgeProb*erandn
    eprobs=0
    do i=1,4
        do j=1,Edges(i)%enum
            eprobs=eprobs+Edges(i)%E(j,5)
            if(eprobs>eprobr) then
                id=Edges(i)%E(j,3)  ! occ
                call det_core(i,j,id)
                if(id>0) then
                    x=Edges(i)%E(j,1)  ! occ
                    y=Edges(i)%E(j,2)  ! occ
                    call get_neighb(x,y)
                    do k=1,3   ! add the detached species onto surface
                        if(Grid(Coor%NB(k,1),Coor%NB(k,2))%occ==0) then
                            call add_species(id,Coor%NB(k,1),Coor%NB(k,2))
                            exit
                        end if
                    end do
                else 
                    call add_species(id)
                end if
                if(mattflag/=0) then  ! update ProbAtt
                    do k=1,MSpecNum
                        MSpecProb=MSpecProb-MSpecList(k)%ProbAtt
                        call mspecies_att_prob(k,MSpecList(k)%ProbAtt)
                        MSpecProb=MSpecProb+MSpecList(k)%ProbAtt
                    end do
                end if
                return
            end if
        end do
    end do
    write(*,*) "error in detach "
  end subroutine 
  !---------------------------------------------------------------------!

  !---------------------------------------------------------------------!

  subroutine update_around(x,y)
    implicit none
    integer,            intent(in) :: x,y   !coordination of this species
    integer                        :: i
    call get_neighb(x,y)
    do i=1, Coor%NBNum
        if(Grid(Coor%NB(i,1),Coor%NB(i,2))%occ>0) then
            call ispecies_prob_update(Grid(Coor%NB(i,1),Coor%NB(i,2))%spec)
        end if
    end do
  end subroutine update_around

  !---------------------------------------------------------------------!

  subroutine species_step()
    implicit none
    double precision  :: probr, probs,randn
    integer           :: id,i
    type(SpecI), pointer  :: p
    call random_number(randn)
    probr=SpecProb*randn
    probs=0
    do id=1,SpecNum                            ! which type
        probs=probs+SpecList(id)%ProbSum
        !write(*,*) "probs", probs
        !write(*,*) "probr", probr
        if(probs>probr) then
            !write(*,*) "species_step"
            call random_number(randn)
            probr=SpecList(id)%ProbSum*randn
            probs=0
            p=>SpecList(id)%head
            do i=1,SpecList(id)%Num            ! which one
                probs=probs+p%ProbSum
                !write(*,*) "p%ProbSum"
                !write(*,*) p%ProbSum
                !write(*,*) "probs", probs
                !write(*,*) "probr", probr
                if(probs>probr) then
                    !write(*,*) "species/id",id,i
                    !write(*,*) "ispecies_step"
                    call ispecies_step(p)
                    p=>null()
                    return
                end if
                if(i/=SpecList(id)%Num) then
                    p=>p%next
                else 
                    call print_err("Error in choosing species step!")
                end if
            end do
        end if
    end do
  end subroutine species_step

  !---------------------------------------------------------------------!

  subroutine ispecies_step(p)
    implicit none
    type(SpecI), pointer, intent(in) :: p
    double precision                 :: probr,probs,randn
    integer                          :: i,j,drct
    integer                          :: id,x,y
    call random_number(randn)
    probr=p%ProbSum*randn
    probs=0
    !if(p%SID==7)   then
    !    write(*,*) "C2H2 ProbSum: ",p%ProbSum
    !    write(*,*) "Desorp: ",p%ProbDesorp
    !    write(*,*) "Decomp: ",p%ProbDecomp
    !    pause
    !end if
    do drct=1, Coor%NBNum           ! drct events which can be diffusion/merge/attach
        if(p%ProbNB(drct)<=0) cycle
        probs=probs+p%ProbNB(drct)
        if(probs>probr) then
            call get_nb_direct(p%Coor(1),p%Coor(2),drct)
            if(Grid(Coor%Drct(1),Coor%Drct(2))%occ==0) then      ! diffusion
                !write(*,*) "diffusion"
                id=p%SID
                x=p%Coor(1)
                y=p%Coor(2)
                call diffusion(id,x,y,drct)
                return
            else if(Grid(Coor%Drct(1),Coor%Drct(2))%occ>0) then  ! merge with species
                !write(*,*) "merge with species"
                !write(*,*) p%SID,p%Coor(1),p%Coor(2),drct
                !write(*,*) "Coor%Drct(1),Coor%Drct(2),Grid"
                !write(*,*) Coor%Drct(1),Coor%Drct(2),Grid(Coor%Drct(1),Coor%Drct(2))%occ
                !write(*,*) "p%Coor,Grid",p%Coor(1),p%Coor(2),Grid(p%Coor(1),p%Coor(2))%occ
                if(Grid(Coor%Drct(1),Coor%Drct(2))%occ==0) call print_err("Grid==0")
                id=p%SID
                x=p%Coor(1)
                y=p%Coor(2)
                call combine(id,x,y,drct=drct)
                return
            else if(Grid(Coor%Drct(1),Coor%Drct(2))%occ==-1) then  ! attach
                id=p%SID
                x=p%Coor(1)
                y=p%Coor(2)
                call att_core(id,Coor%Drct(1),Coor%Drct(2))  ! first att, otherwise Coor%Drct might change
                call del_species(id,x,y)
                  
                if(mattflag/=0) then  ! update ProbAtt due to the att_core
                    do j=1,MSpecNum
                        MSpecProb=MSpecProb-MSpecList(j)%ProbAtt
                        call mspecies_att_prob(j,MSpecList(j)%ProbAtt)
                        MSpecProb=MSpecProb+MSpecList(j)%ProbAtt
                    end do
                end if

                return
            !else occ==-10, do nothing
            end if
        end if
     end do

     do i=1, SpecEvt(p%SID)%MergeNum                  ! merge with mspecies
        if(p%ProbMerge(i)<=0) cycle
        probs=probs+p%ProbMerge(i)
        if(probs>probr) then
            !write(*,*) "merge with mspecies"
            id=p%SID
            x=p%Coor(1)
            y=p%Coor(2)
            call combine(id,x,y,evt=i)
            return
        end if
     end do

     do i=1, SpecEvt(p%SID)%DecompNum                  ! decomposition
        if(p%ProbDecomp(i)<=0) cycle
        probs=probs+p%ProbDecomp(i)
        if(probs>probr) then
            !write(*,*) "decomp"
            id=p%SID
            x=p%Coor(1)
            y=p%Coor(2)
            call decomp(id,i,x,y)
            return
        end if
     end do

     if(p%ProbDesorp>0) then                   ! desorption
         probs=probs+p%ProbDesorp
         if(probs>probr) then
             !write(*,*) "desorp"
             id=p%SID
             x=p%Coor(1)
             y=p%Coor(2)
             call del_species(id,x,y)
             call update_event_num("desorption",DesorpI(id))
             return
         end if
     end if

     call species_prob_update()
     !write(*,*) "p: ",p%SID
     !write(*,*) "probs,Probsum,probNB: ",probs,p%ProbSum,p%probNB
     !if(SpecEvt(i)%MergeNum>0)  write(*,*) "merge: ",p%ProbMerge 
     !if(SpecEvt(i)%DecompNum>0) write(*,*) "Decomp: ",p%ProbDecomp
     !if(p%ProbDesorp>0)         write(*,*) "Desorp: ",p%ProbDesorp
     !call print_err("Error in ispecies_step")

  end subroutine ispecies_step

  !---------------------------------------------------------------------!

  subroutine species_prob_update()
    implicit none
    type(SpecI), pointer       :: p
    integer        :: i,num
    SpecProb=0
    do i=1,SpecNum
        SpecList(i)%ProbSum=0
        if(associated(SpecList(i)%head)) then
            num=0
            p=>SpecList(i)%head
            do
                call iispecies_prob_update(p)
                SpecList(i)%ProbSum=SpecList(i)%ProbSum+p%ProbSum
                num=num+1
                if(.not.associated(p%next)) exit;
                p=>p%next
            end do
            if(num/=SpecList(i)%Num) call print_err("Error in species_prob_update",3)
        end if
        SpecProb=SpecProb+SpecList(i)%ProbSum
     end do
     p=>null()
  end subroutine species_prob_update

  !---------------------------------------------------------------------!

  subroutine ispecies_prob_update(p)
    implicit none
    type(SpecI), pointer, intent(in)  :: p
    if(p%ProbSum>0) then
        SpecProb=SpecProb-p%ProbSum
        SpecList(p%SID)%ProbSum=SpecList(p%SID)%ProbSum-p%ProbSum
    end if
    call iispecies_prob_update(p)
    SpecProb=SpecProb+p%ProbSum
    SpecList(p%SID)%ProbSum=SpecList(p%SID)%ProbSum+p%ProbSum
  end subroutine ispecies_prob_update

  !---------------------------------------------------------------------!

  subroutine iispecies_prob_update(p)
    implicit none
    type(SpecI), pointer,  intent(in)  :: p
    integer           :: voidnum
    integer           :: i,drct
    integer           :: id1, id2
    integer           :: nn !for attach, num of neighboring sites
    voidnum=0
    p%ProbSum=0
    do drct=1, Coor%NBNum
        call get_nb_direct(p%Coor(1),p%Coor(2),drct=drct)
        !write(*,*) "Coor%Drct"
        !write(*,*) Coor%Drct
        if(Grid(Coor%Drct(1),Coor%Drct(2))%occ==0) then  ! diffusion
            voidnum=voidnum+1
            if(SpecEvt(p%SID)%DiffProb>0) then
                p%ProbNB(drct)=SpecEvt(p%SID)%DiffProb
                p%ProbSum=p%ProbSum+p%ProbNB(drct)
            else
                p%ProbNB(drct)=0
            end if
        else if(Grid(Coor%Drct(1),Coor%Drct(2))%occ>0) then ! merge with species
            id1=p%SID
            id2=Grid(Coor%Drct(1),Coor%Drct(2))%occ
            if(MergeMProb(id1,id2)>0) then
                p%ProbNB(drct)=MergeMProb(id1,id2)
                p%ProbSum=p%ProbSum+p%ProbNB(drct)
            else
                p%ProbNB(drct)=0
            end if
        else if(Grid(Coor%Drct(1),Coor%Drct(2))%occ==-1) then ! attach
            call get_nn(Coor%Drct(1),Coor%Drct(2),nn)
            if(SpecEvt(p%SID)%AttProb(nn)>0) then
                p%ProbNB(drct)=SpecEvt(p%SID)%AttProb(nn)
                p%ProbSum=p%ProbSum+p%ProbNB(drct)
            else
                p%ProbNB(drct)=0
            end if
        else
            p%ProbNB(drct)=0    ! no event, no summation
        end if
    end do

    if(voidnum>0) then  !decomp & merge with mspecies
        do i=1, SpecEvt(p%SID)%MergeNum                    ! merge with mspecies
            id1=SpecEvt(p%SID)%MergeRI(i,1)
            !p%ProbMerge(i)=SpecEvt(p%SID)%MergeProbP(i)*voidnum*MSpecList(-id1)%Num/SiteNum
            p%ProbMerge(i)=SpecEvt(p%SID)%MergeProbP(i)*MSpecList(-id1)%Num/SiteNum
            p%ProbSum=p%ProbSum+p%ProbMerge(i)
        end do
        if(allocated(p%ProbDecomp)) then                    ! decomp
            p%ProbDecomp(:)=SpecEvt(p%SID)%DecompProbP(:)
            do i=1,SpecEvt(p%SID)%DecompNum             
                p%ProbSum=p%ProbSum+p%ProbDecomp(i)
            end do
        end if
    else 
        if(allocated(p%ProbMerge))  p%ProbMerge=0
        if(allocated(p%ProbDecomp)) p%ProbDecomp=0
    end if

    ! p%ProbDesorp do not need to update
    if(p%ProbDesorp>0) p%ProbSum=p%ProbSum+p%ProbDesorp    !desorp

    !write(*,*) "id,prob: ",p%SID,p%ProbSum,p%ProbNB,p%ProbDecomp

  end subroutine iispecies_prob_update

  !---------------------------------------------------------------------!

  subroutine mspecies_step()
    implicit none
    double precision  :: probr, probs,randn
    integer           :: i,j,k,id,nn,m
    call random_number(randn)
    probr=MSpecProb*randn
    probs=0
    do i=1,MSpecNum
        if(MSpecList(i)%ProbDesorp>0) then
            probs=probs+MSpecList(i)%ProbDesorp
            if(probs>probr) then
                call del_species(MSpecList(i)%ID)
                call update_event_num("desorption",DesorpI(SpecNum-MSpecList(i)%ID))
                return
            end if
        end if
        do j=1,MSpecEvt(i)%DecompNum
            if(MSpecList(i)%ProbDecomp(j)>0) then
                probs=probs+MSpecList(i)%ProbDecomp(j)
                if(probs>probr) then
                    call decomp(MSpecList(i)%ID,j)
                    return
                end if
            end if
        end do
        if(mattflag/=0) then 
            do j=1,4
               do k=1,Edges(j)%enum
                  if(Edges(j)%E(k,3)==0) then
                      nn=Edges(j)%E(k,4)
                      probs=probs+MSpecEvt(i)%AttProb(nn)*MSpecList(i)%Num/SiteNum
                      if(probs>probr) then
                          call att_core(-i,Edges(j)%E(k,1),Edges(j)%E(k,2))  ! first att, otherwise Coor%Drct might change
                          if(mattflag/=0) then  ! update ProbAtt due to the att_core
                              do m=1,MSpecNum
                                  MSpecProb=MSpecProb-MSpecList(m)%ProbAtt
                                  call mspecies_att_prob(j,MSpecList(m)%ProbAtt)
                                  MSpecProb=MSpecProb+MSpecList(m)%ProbAtt
                              end do
                          end if
                          call del_species(-i)
                          return
                      end if
                  end if
               end do
            end do
        end if
    end do
    call print_err("Error in choosing mspecies step!")
  end subroutine mspecies_step

  !---------------------------------------------------------------------!

  subroutine mspecies_att_prob(id,probatt)
      implicit none
      integer, intent(in)  :: id
      double precision, intent(out) :: probatt
      integer :: i,j,nn
      probatt=0
      do i=1,4
         do j=1,Edges(i)%enum
            if(Edges(i)%E(j,3)==0) then
                nn=Edges(i)%E(j,4)
                probatt=probatt+MSpecEvt(id)%AttProb(nn)
            end if
         end do
      end do
      probatt=probatt*MSpecList(id)%Num/SiteNum
  end subroutine 

  !---------------------------------------------------------------------!

  subroutine mspecies_prob_update(id)
    implicit none
    integer, optional, intent(in) :: id
    integer                       :: i,j,nn
    double precision              :: prob
    
    if(present(id)) then
        if(MSpecEvt(-id)%DesorpBrr>0) then
            MSpecProb=MSpecProb-MSpecList(-id)%ProbDesorp
            MSpecList(-id)%ProbDesorp=MSpecEvt(-id)%DesorpProbP*MSpecList(-id)%Num
            MSpecProb=MSpecProb+MSpecList(-id)%ProbDesorp
        end if
        if(MSpecEvt(-id)%DecompNum>0) then
            MSpecProb=MSpecProb-sum(MSpecList(-id)%ProbDecomp)
            MSpecList(-id)%ProbDecomp(:)=MSpecEvt(-id)%DecompProbP(:)*MSpecList(-id)%Num
            MSpecProb=MSpecProb+sum(MSpecList(-id)%ProbDecomp)
        end if
        if(mattflag/=0) then   ! mspecies attach onto core
            MSpecProb=MSpecProb-MSpecList(-id)%ProbAtt
            call mspecies_att_prob(-id,prob)
            MSpecList(-id)%ProbAtt=prob
            MSpecProb=MSpecProb+MSpecList(-id)%ProbAtt
        end if
        return
    end if

    MSpecProb=0
    do i=1, MSpecNum
        if(MSpecEvt(i)%DesorpBrr>0) then
            MSpecList(i)%ProbDesorp=MSpecEvt(i)%DesorpProbP*MSpecList(i)%Num
            MSpecProb=MSpecProb+MSpecList(i)%ProbDesorp
        end if
        if(MSpecEvt(i)%DecompNum>0) then
            if(.not.allocated(MSpecList(i)%ProbDecomp)) then
                allocate(MSpecList(i)%ProbDecomp(MSpecEvt(i)%DecompNum))
            end if
            MSpecList(i)%ProbDecomp(:)=MSpecEvt(i)%DecompProbP(:)*MSpecList(i)%Num
            MSpecProb=MSpecProb+sum(MSpecList(i)%ProbDecomp)
        end if
        if(mattflag/=0) then ! mspecies attach onto core
            MSpecProb=MSpecProb-MSpecList(i)%ProbAtt
            call mspecies_att_prob(i,prob)
            MSpecList(i)%ProbAtt=prob
            MSpecProb=MSpecProb+MSpecList(i)%ProbAtt
        end if
    end do

  end subroutine mspecies_prob_update

  !---------------------------------------------------------------------!

  subroutine mmmerge_step()
    implicit none
    double precision  :: probr, probs,randn
    integer           :: i,id
    call random_number(randn)
    probr=MMMergeProb*randn
    probs=0
    do i=1,MMergeNum
        probs=probs+MMergeProb(i)
        if(probs>probr) then
            call del_species(MMergeRI(i,1))
            call del_species(MMergeRI(i,2))
            call add_species(MMergeRI(i,3))
            call update_event_num("merge",MMergeRI(i,4))
            return
        end if
    end do
    call print_err("Error in choosing mspecies step!")
  end subroutine mmmerge_step

  !---------------------------------------------------------------------!

  subroutine mmmerge_prob_update()
    implicit none
    integer           :: id1,id2
    integer           :: i
    MMMergeProb=0
    do i=1, MMergeNum
        id1=MMergeRI(i,1)
        id2=MMergeRI(i,2)
        MMergeProb(i)=MMergeProbP(i)*MSpecList(-id1)%Num*MSpecList(-id2)%Num/SiteNum
        MMMergeProb=MMMergeProb+MMergeProb(i)
    end do
  end subroutine mmmerge_prob_update

  !---------------------------------------------------------------------!

  subroutine print_species_info(u)
    implicit none
    integer,intent(in) :: u
    integer            :: i,j,k

    write(u,*) "AllSpecNum: ",io_adjustl(AllSpecNum)
    write(u,*) "SpecNum:    ",io_adjustl(SpecNum)
    if(SpecNum>0) then
        write(u,*) "SpecName:   ",SpecName
        write(u,'(" ID: "$)')
        do i=1, SpecNum
            write(u,'(I10$)') i
        end do
        write(u,*)
    end if

    write(u,*) "MSpecNum:   ",io_adjustl(MSpecNum)
    if(MSpecNum>0) then
        write(u,*) "MSpecName:  ",MSpecName
        write(u,'(" ID: "$)')
        do i=1, MSpecNum
            write(u,'(I10$)') -i
        end do
        write(u,*)
    end if
    write(u,*)

    write(u,*) "MergeM"
    do i=1, AllSpecNum
        do j=1, AllSpecNum
            if(MergeMBrr(i,j)>0) then
                write(u,'(" ",I2$)') MergeM(i,j)
            else
                write(u,'(A$)') " --"
            end if
        end do
        write(u,*)
    end do
    write(u,*)

    write(u,*) "MergeMI"
    do i=1, AllSpecNum
        do j=1, AllSpecNum
            if(MergeMBrr(i,j)>0) then
                write(u,'(" ",I2$)') MergeMI(i,j)
            else
                write(u,'(A$)') " --"
            end if
        end do
        write(u,*)
    end do
    write(u,*)

    write(u,*) "MergeMBrr"
    do i=1, AllSpecNum
        do j=1, AllSpecNum
            if(MergeMBrr(i,j)>0) then
                write(u,'(" ",F5.3$)') MergeMBrr(i,j)
            else
                write(u,'(A$)') "  --- "
            end if
        end do
        write(u,*)
    end do
    write(u,*)
    write(u,*) "MergeMProb"
    do i=1, AllSpecNum
        do j=1, AllSpecNum
            if(MergeMProb(i,j)>0) then
                write(u,'(" ",E9.2$)') MergeMProb(i,j)
            else
                write(u,'(A$)') "   ----   "
            end if
        end do
        write(u,*)
    end do
    write(u,*)

    if(MMergeNum>0) then
        write(u,'(A,I3)') " MMergeNum:  ",MMergeNum
        do i=1, MMergeNum
            write(u,'(A,I3,A)') " MMerge Reaction Index ", MMergeRI(i,4),&
                              "                      -----------------------------"
            write(u,'(A,3(I3," "))') " MMergeRI:   ",MMergeRI(i,1:3)
            write(u,'(A,F6.3)') " MMergeBrr:  ",MMergeBrr(i)
            write(u,'(A,E9.2)') " MMergeProb: ",MMergeProbP(i)
        end do
    end if
    write(u,*)

    if(SpecNum>0) then
        write(u,*) repeat('=',70)
        write(u,*) io_center("SpecE", 70)
        write(u,*) repeat('=',70)
    end if
    do i=1, SpecNum
        write(u,*) "DiffBrr/DiffProb/DesorpBrr/DesorpProb/AttBrr/AttProb/DetBrr/DetProb/DecompNum"
        if(SpecEvt(i)%DiffBrr>0) then
            write(u,'(F6.3," ",E9.2," "$)') SpecEvt(i)%DiffBrr,SpecEvt(i)%DiffProb
        else 
            write(u,'(A$)') "  ---    -----   "
        end if
        if(SpecEvt(i)%DesorpBrr>0) then
            write(u,'(F6.3," ",E9.2," "$)') SpecEvt(i)%DesorpBrr,SpecEvt(i)%DesorpProbP
        else 
            write(u,'(A$)') "  ---    -----   "
        end if
        if(SpecEvt(i)%AttBrr(1)>0) then
            write(u,'(F6.3," ",E9.2," "$)') SpecEvt(i)%AttBrr(1),SpecEvt(i)%AttProb(1)
        else 
            write(u,'(A$)') "  ---    -----   "
        end if
        if(SpecEvt(i)%DetBrr(1)>0) then
            write(u,'(F6.3," ",E9.2," "$)') SpecEvt(i)%DetBrr(1),SpecEvt(i)%DetProb(1)
        else 
            write(u,'(A$)') "  ---    -----   "
        end if
        write(u,'(I3)') SpecEvt(i)%DecompNum
        do j=1,SpecEvt(i)%DecompNum
            write(u,'(A,I3,A)') " Decomp Reaction Index ", SpecEvt(i)%DecompRI(j,3),&
                              "                      -----------------------------"
            write(u,'(A,3(I3," "))') " DecompRI", SpecEvt(i)%DecompRI(j,1:2)
            write(u,'(A,F6.3)') " DecompBrr",     SpecEvt(i)%DecompBrr(j)
            write(u,'(A,E9.2)') " DecompProb",    SpecEvt(i)%DecompProbP(j)
        end do
        do j=1,SpecEvt(i)%MergeNum
            write(u,'(A,I3,A)') " Merge Reaction Index (with MSpecies only) ", SpecEvt(i)%MergeRI(j,3),&
                              "  -----------------------------"
            write(u,'(A,3(I3," "))') " MergeRI", SpecEvt(i)%MergeRI(j,1:2)
            write(u,'(A,F6.3)') " MergeBrr",     SpecEvt(i)%MergeBrr(j)
            write(u,'(A,E9.2)') " MergeProb",    SpecEvt(i)%MergeProbP(j)
        end do
        write(u,*)
    end do

    if(MSpecNum>0) then
        write(u,*) repeat('=',70)
        write(u,*) io_center("MSpecE", 70)
        write(u,*) repeat('=',70)
    end if
    do i=1, MSpecNum
        write(u,*) "DesorpBrr/DesorpProb/DecompNum"
        if(MSpecEvt(i)%DesorpBrr>0) then
            write(u,'(F6.3,"   ",E9.2,"   "$)') MSpecEvt(i)%DesorpBrr,MSpecEvt(i)%DesorpProbP
        else 
            write(u,'(A$)') "  ---      -----     "
        end if
        write(u,'(I3)') MSpecEvt(i)%DecompNum
        do j=1,MSpecEvt(i)%DecompNum
            write(u,'(A,I3,A)') " Decomp Reaction Index ", MSpecEvt(i)%DecompRI(j,3),&
                              "                      -----------------------------"
            write(u,'(A,3(I3," "))') " DecompRI", MSpecEvt(i)%DecompRI(j,1:2)
            write(u,'(A,F6.3)') " DecompBrr",     MSpecEvt(i)%DecompBrr(j)
            write(u,'(A,E9.2)') " DecompProb",    MSpecEvt(i)%DecompProbP(j)
        end do
        write(u,*)
    end do

  end subroutine print_species_info

  !---------------------------------------------------------------------!

  subroutine save_species(u)
    implicit none
    integer, intent(in)          :: u
    integer                      :: i
    if(SpecNum>0) then
        write(u,'(2A)') "NumberOfSpeciesType  ", trim(io_adjustl(SpecNum))
        write(u,'(A)') "NumberOfEachSpecies "
        do i=1, SpecNum
            write(u,'(A," "$)') trim(io_adjustl(SpecList(i)%Num))
        end do
        write(u,*)
        write(u,*)
    end if
    if(MSpecNum>0) then
        write(u,'(2A)') "NumberofMSpeciesType ", trim(io_adjustl(MSpecNum))
        write(u,'(A)') "NumberofEachMSpecies "
        do i=1, MSpecNum
            write(u,'(A," "$)') trim(io_adjustl(MSpecList(i)%Num))
        end do
        write(u,*)
        write(u,*)
    end if

  end subroutine save_species

  !---------------------------------------------------------------------!
  subroutine print_species_name(u)
    implicit none
    integer, intent(in)          :: u
    integer                      :: i
    do i=1, MSpecNum
        write(u,"(' ',A8,'  '$)") trim(adjustl(MSpecName(i)))
    end do
    do i=1, SpecNum
        write(u,"(' ',A8,'  '$)") trim(adjustl(SpecName(i)))
    end do
    write(u,*)
  end subroutine print_species_name

  !---------------------------------------------------------------------!


  subroutine print_species_density(u,num)
    implicit none
    integer, intent(in)          :: u
    integer,optional,intent(in)  :: num
    integer                      :: i
    do i=1, MSpecNum
        if(present(num)) then
            write(u,'(I10," "$)')   MSpecList(i)%Num
        else
            write(u,'(" ",E10.3$)') real(MSpecList(i)%Num)/SiteNum
        end if
    end do
    do i=1, SpecNum
        if(present(num)) then
            write(u,'(I10," "$)')   SpecList(i)%Num
        else
            write(u,'(" ",E10.3$)') real(SpecList(i)%Num)/SiteNum
        end if
    end do
    write(u,*)
  end subroutine print_species_density

end module species
