!-----------------------------------------------------------------------
!                   input.f90 - input file handling
!-----------------------------------------------------------------------
! This file is part of the TEAKS package.
!
! Copyright (C) 2020-2026 Pai Li
! 2020-01-27 Pai Li
!-----------------------------------------------------------------------

module input

  use teio,     only: teio_readline,           &
                      teio_assert_file_exists, &
                      PATHLEN, LINELEN, SPCSLEN

  use io,       only: io_adjustl,  &
                      io_lower,    &
                      io_readnext, &
                      io_readval,  &
                      io_unit

  implicit none
  private
  save

  public  :: read_InpTE,        &
             read_InpKMC,       &
             inp_print_info,    &
             print_err

  private :: inp_read_value,    &
             inp_read_value_c1, &
             inp_read_value_d1, &
             inp_read_value_i1, &
             inp_read_value_i2, &
             inp_read_value_l1, &
             inp_read_array,    &
             inp_read_array_c,  &
             inp_read_array_i,  &
             inp_read_array_d,  &
             inp_read_block_im, & 
             inp_read_block_cm, & 
             inp_read_block_ca, & 
             inp_read_block,    &
             read_long_int


  !--------------------------------------------------------------------!
  !                      generic input data type                       !
  !--------------------------------------------------------------------!

  type, public :: InputData

     !-----------------------------------------------------------------!
     ! init             .true. when memory has been allocated          !
     ! file             name of the input file the data was read from  !
     ! outFileName      name of the output file                        !
     ! mode             run mode (predict, mc)                         !
     ! verbosity        verbosity level (0=low, 1=normal, 2=high)      !
     !                                                                 !
     ! typeName(i)      name of i-th atomic species                    !
     ! atomicEnergy(i)  atomic energy of i-th atomic species           !
     !                                                                 !
     ! activeType(i)    == 1, if type i is "active"                    !
     !                                                                 !
     ! netFile(i)       path to NN potential file for species i        !
     ! netArch(i)       architecture (as string) of NN for species i   !
     !                                                                 !
     ! do_forces        .true. if forces shall be calculated           !
     ! do_timing        .true. if timings shall be saved               !
     ! print_atomic_energies  if .true., atomic energies will be saved !
     !                                                                 !
     ! nStrucs          number of structures to run calculations for   !
     ! strucFile(i)     path to i-th atomic structure file             !
     !                                                                 !
     !-------------------------- Monte-Carlo --------------------------!
     ! T                Temperature                                    !
     ! VibMode          the way to calculate vibration frequency.      !
     !                  1 for KT/h,                                    !
     ! Vib              Vibration frequency                            !
     ! MSize            Mesh size of the substrate                     !
     ! CoreSize         Core size of the crystal island                !
     ! MType            Mesh type. 1 for square, 2 for hexagon grid    !
     ! TStep            Total Steps                                    !
     ! PrtDStep         Interval Steps for Density Print               !
     ! PrtEStep         Interval Steps for Event Count Print           !
     !                                                                 !
     ! GasNum           Number of gas' types                           !
     ! GasName(i)        i-th Gas' name                                !
     ! GasPrs(i)        Gas Pressure                                   !
     ! GasAdsorbRate(i) Gas adsorption rate                            !
     ! GasProd(i,2)     adsorption product                             !
     !                                                                 !
     ! SpecNum          Number of species types                        !
     ! SpecName(i)       i-th Species' name                            !
     ! IniSpecDens(i)   Initial species densities                      !
     !                                                                 !
     ! DiffBrr(i)       Diffusion barrier                              !
     ! AttBrr(i)        Attachment barrier                             !
     ! DetBrr(i)        Detachment barrier                             !
     !                                                                 !
     ! MergeNum         Merge event number                             !
     ! MergeR(i,3)      Merge event: R1, R2, P                         !
     ! MergeBrr(i)      Merge event barrier                            !
     !                                                                 !
     ! DecompNum        Decomp event number                            !
     ! DecompR(i,3)     Decomp event: R, P1, P2                        !
     ! DecompBrr(i)     Decomp event barrier                           !
     !                                                                 !
     ! DesorpNum        Desorp event number                            !
     ! DesorpR(i)       Desorp event: R                                !
     ! DesorpBrr(i)     Desorp event barrier                           !
     !                                                                 !
     ! MSpecName         Species' name with MFA                        !
     !                                                                 !
     !-----------------------------------------------------------------!

     logical                                             :: init = .false.
     character(len=PATHLEN)                              :: file
     character(len=SPCSLEN)                              :: Mode
     character(len=5)                                    :: Restart

     ! KMC simulation parameters
     double precision                                    :: T   !temperature
     integer                                             :: VibMode 
     double precision                                    :: Vib
     integer                                             :: MeshType
     integer,                dimension(2)                :: MeshSize
     logical                                             :: FlagCore
     integer,                dimension(2)                :: CoreSize
     integer                                             :: FlagFixDen
     integer,                dimension(2)                :: FixDenArea
     integer                                             :: MisMatchN
     integer,                dimension(:),   allocatable :: MisMatchP


     double precision                                    :: StepRead
     integer(kind=8)                                     :: TotStep
     integer(kind=8)                                     :: PrtSpecDenStep!Step inteval for species-info print
     integer(kind=8)                                     :: SaveStatusStep

     integer                                             :: GasNum
     character(len=SPCSLEN), dimension(:),   allocatable :: GasName
     double precision,       dimension(:),   allocatable :: GasPrs
     double precision,       dimension(:),   allocatable :: GasAdsorbRate
     character(len=SPCSLEN), dimension(:,:), allocatable :: GasProd
     
     integer                                             :: MSpecNum
     character(len=SPCSLEN), dimension(:),   allocatable :: MSpecName
     double precision,       dimension(:),   allocatable :: MIniSpecDen
     integer,                dimension(:),   allocatable :: MIniSpecNum
     integer,                dimension(:),   allocatable :: MSpecDenFix

     integer                                             :: mattflag
     double precision,       dimension(:,:), allocatable :: MAttBrr
     double precision,       dimension(:,:), allocatable :: MDetBrr

     integer                                             :: SpecNum
     character(len=SPCSLEN), dimension(:),   allocatable :: SpecName
     double precision,       dimension(:),   allocatable :: IniSpecDen
     integer,                dimension(:),   allocatable :: IniSpecNum
     integer,                dimension(:),   allocatable :: SpecDenFix
     ! index of species is positive integer

     integer                                             :: attflag
     double precision,       dimension(:),   allocatable :: DiffBrr
     double precision,       dimension(:,:), allocatable :: AttBrr
     double precision,       dimension(:,:), allocatable :: DetBrr

     integer                                             :: MergeNum
     character(len=SPCSLEN), dimension(:,:), allocatable :: MergeR  !MergeReaction: R1, R2, P
     double precision,       dimension(:),   allocatable :: MergeBrr

     integer                                             :: DecompNum
     character(len=SPCSLEN), dimension(:,:), allocatable :: DecompR !DecompReaction: R1, P1, P2
     double precision,       dimension(:),   allocatable :: DecompBrr

     integer                                             :: DesorpNum
     character(len=SPCSLEN), dimension(:),   allocatable :: DesorpR !Desorption: R
     double precision,       dimension(:),   allocatable :: DesorpBrr

     ! for restart
     integer                                             :: GridOccNum
     integer,                dimension(:,:), allocatable :: GridOcc
     integer                                             :: SpecNumR 
     integer,                dimension(:),   allocatable :: NumSpecies
     integer                                             :: MSpecNumR
     integer,                dimension(:),   allocatable :: NumMSpecies
     integer,                dimension(4)                :: EdgeSiteNum
     integer,                dimension(:),   allocatable :: Edge1
     integer,                dimension(:),   allocatable :: Edge2
     integer,                dimension(:),   allocatable :: Edge3
     integer,                dimension(:),   allocatable :: Edge4
     
  end type InputData

  !--------------------------------------------------------------------!

  interface inp_read_value
     module procedure inp_read_value_c1, inp_read_value_d1, &
                      inp_read_value_i1, inp_read_value_i2, inp_read_value_l1
  end interface inp_read_value

  interface inp_read_array
     module procedure inp_read_array_i, inp_read_array_c, &
                      inp_read_array_d
  end interface inp_read_array

  interface inp_read_block
     module procedure inp_read_block_im, inp_read_block_cm, &
                      inp_read_block_ca
  end interface inp_read_block

  !----------------------- return status values -----------------------!

  integer, parameter, private :: S_OK    = 0
  integer, parameter, private :: S_ERROR = 2
  integer, parameter, private :: S_NOT   = 3

  double precision,parameter,public :: Kb=8.617333262145E-5 !Boltzmann constant in eV*K-1
  double precision,parameter,public :: h=4.135667696E-15    !Plank constant in eV*s
contains

  !--------------------------------------------------------------------!
  !                         default parameters                         !
  !--------------------------------------------------------------------!

  subroutine inp_defaults(inp)

    implicit none

    type(InputData), intent(inout) :: inp

    inp%VibMode=1
    inp%Restart="False"
    inp%TotStep=1000
    inp%PrtSpecDenStep=100
    inp%SaveStatusStep=100

    inp%GasNum=0
    inp%SpecNum=0
    inp%MSpecNum=0
    inp%MergeNum=0
    inp%DecompNum=0
    inp%DesorpNum=0


  end subroutine inp_defaults

  !====================================================================!
  !                                                                    !
  !                       parsers for input file                       !
  !                                                                    !
  !====================================================================!

  !------------------------------- KMC --------------------------------!

  function read_InpKMC(file) result(inp)

    implicit none

    character(len=*), intent(in) :: file
    type(InputData)              :: inp

    integer                :: iline, stat
    integer                :: u
    integer                :: i,j
    logical                :: fexists

    call teio_assert_file_exists(file)
    inp%file = trim(file)
    call inp_defaults(inp)

    u = io_unit()
    open(u, file=trim(adjustl(file)), status='old', action='read')

    call inp_read_value(u, 'mode',         inp%Mode)
    call inp_read_value(u, 'restart',      inp%Restart)
    call inp_read_value(u, 'temperature',  inp%T)
    call inp_read_value(u, 'vibmode',      inp%VibMode)

    if(inp%VibMode==1)  then
      inp%Vib=Kb*inp%T/h
    else
      call inp_read_value(u, 'vib',          inp%Vib,stat)
      if(stat>0)  call print_err("Error in reading inp%Vib")
    end if
    call inp_read_value(u, 'meshtype',     inp%MeshType)
    call inp_read_array(u, 'meshsize',     inp%MeshSize)
    call inp_read_array(u, 'coresize',     inp%CoreSize, stat)
    inp%FlagCore=.False.
    if ( stat==0 ) then
        inp%FlagCore=.True.
        if(inp%CoreSize(1)>inp%MeshSize(1)*0.9.or.inp%CoreSize(2)>inp%MeshSize(2)*0.9) then
            call print_err("CoreSize should be smaller than 0.9*MeshSize")
        end if
        if(inp%MeshType==1) then
            if(mod(inp%CoreSize(1),2)/=0 .or. mod(inp%CoreSize(2),2)/=0) then
                call print_err("CoreSize should be even numbers for Hexagonal grid")
            end if
        end if
    end if
    call inp_read_value(u, 'mismatchn',      inp%MisMatchN, stat)
    if(stat==0) then
        if(inp%MisMatchN<1) call print_err("MisMatch should be larger than 0")
        allocate(inp%MisMatchP(inp%MisMatchN))
        call inp_read_array(u, 'mismatchp',     inp%MisMatchP)
    end if
        

    call inp_read_array(u, 'fixdenarea',     inp%FixDenArea, stat)
    inp%FlagFixDen=0
    if ( stat==0 ) then
        inp%FlagFixDen=.True.
        if(inp%FixDenArea(1)<=inp%CoreSize(1)) then
            call print_err("error: inp%FixDenArea(1)<=inp%CoreSize(1))")
        end if
        if(inp%FIxDenArea(2)/=inp%MeshSize(2)) then
            call print_err("error: inp%FIxDenArea(2)\=inp%MeshSize(2)")
        end if
    end if

    call inp_read_value(u, 'totalstep',    inp%StepRead,stat)
    if(stat==0) inp%TotStep=read_long_int(inp%StepRead)
    call inp_read_value(u, 'printspecdenstep', inp%StepRead,stat)
    if(stat==0) then
        inp%PrtSpecDenStep=read_long_int(inp%StepRead)
    else
        inp%PrtSpecDenStep=inp%TotStep/10
    end if
    call inp_read_value(u, 'savestatus',       inp%StepRead,stat)
    if(stat==0) then
        inp%SaveStatusStep=read_long_int(inp%StepRead)
    else
        inp%SaveStatusStep=inp%TotStep/10
    end if

  !-------------------------------Gas----------------------------------!
    call inp_read_value(u, 'gasnumber',inp%GasNum, stat)
    if ( stat==0 .and. inp%GasNum/=0 ) then
        allocate(inp%GasName(inp%GasNum),inp%GasPrs(inp%GasNum))
        allocate(inp%GasAdsorbRate(inp%GasNum),inp%GasProd(inp%GasNum,3))
        call inp_read_array(u, 'gasname',               inp%GasName, stat)
        if ( stat>0 ) call print_err("Error in reading GasName")
        call inp_read_array(u, 'gaspressure',           inp%GasPrs,  stat)
        if ( stat>0 ) call print_err("Error in reading GasPressure")
        call inp_read_array(u, 'gasadsorptionrate',inp%GasAdsorbRate,  stat)
        if ( stat>0 ) call print_err("Error in reading GasAdsorptionRate")
        call inp_read_block(u, 'gasadsorbreaction', inp%GasProd,stat=stat)
        if ( stat>0 ) call print_err("Error in reading GasAdsorbReaction")
    end if

  !------------------------------MSpecies-------------------------------!

    call inp_read_value(u, 'mspeciesnumber', inp%MSpecNum,  stat)
    if ( stat==0 .and. inp%MSpecNum>0 ) then   !S_OK==0
        allocate(inp%MSpecName(inp%MSpecNum),inp%MIniSpecDen(inp%MSpecNum))
        allocate(inp%MIniSpecNum(inp%MSpecNum))
        allocate(inp%MSpecDenFix(inp%MSpecNum))
        call inp_read_array(u, 'mspeciesname',    inp%MSpecName, stat)
        if ( stat > 0 ) call print_err("Error in reading MSpeciesName")
        inp%MIniSpecDen=0 ! default value for initial density of species
        inp%MIniSpecNum=0 ! default value for initial density of species
        call inp_read_array(u, 'initialmspeciesdensity', inp%MIniSpecDen)
        call inp_read_array(u, 'initialmspeciesnumber',     inp%MIniSpecNum)
        inp%MSpecDenFix=0
        call inp_read_array(u, 'fixmspeciesdensity',   inp%MSpecDenFix)

        ! MSpecies attach and detach
        inp%mattflag=0
        if (inp%FlagCore.eqv..True. ) then
            allocate(inp%MAttBrr(inp%MSpecNum,3),inp%MDetBrr(inp%MSpecNum,3))
            inp%MAttBrr=-1
            inp%MDetBrr=-1
            call inp_read_array(u, 'mattachbarrier',     inp%MAttBrr(:,1), stat)
            if ( stat > 0 ) then
                inp%mattflag=0 !1 for all, 2 for only attachbarrier1 provided, 0 for nothing
            else 
                call inp_read_array(u, 'mattachbarrier2',    inp%MAttBrr(:,2), stat)
                if ( stat > 0 ) then
                    inp%MAttBrr(:,2)=inp%MAttBrr(:,1)
                    inp%MAttBrr(:,3)=inp%MAttBrr(:,1)
                    inp%mattflag=2
                else
                    inp%mattflag=1
                    call inp_read_array(u, 'mattachbarrier3',    inp%MAttBrr(:,3), stat)
                    if ( stat > 0 ) call print_err("mattachbarrier3 should be provided")
                end if 
            end if
            if(inp%mattflag/=0) then
                call inp_read_array(u, 'mdetachbarrier',     inp%MDetBrr(:,1), stat)
                if(inp%mattflag==1) then
                    call inp_read_array(u, 'mdetachbarrier2',    inp%MDetBrr(:,2), stat)
                    if ( stat > 0 ) inp%MDetBrr(:,2)=inp%MDetBrr(:,1)
                    call inp_read_array(u, 'mdetachbarrier3',    inp%MDetBrr(:,3), stat)
                    if ( stat > 0 ) inp%MDetBrr(:,3)=inp%MDetBrr(:,1)
                end if
            end if
        end if
    end if
  !------------------------------Species-------------------------------!
    call inp_read_value(u, 'speciesnumber',         inp%SpecNum, stat)
    if ( stat==0 .and. inp%SpecNum/= 0 ) then   !S_OK==0
        allocate(inp%SpecName(inp%SpecNum))
        allocate(inp%IniSpecDen(inp%SpecNum))
        allocate(inp%IniSpecNum(inp%SpecNum))
        allocate(inp%SpecDenFix(inp%SpecNum))
        allocate(inp%DiffBrr(inp%SpecNum))
        call inp_read_array(u, 'speciesname',           inp%SpecName, stat)
        if ( stat > 0 ) call print_err("Error in reading SpeciesName")

        inp%IniSpecDen=0 ! default value for initial density of species
        inp%IniSpecNum=0 ! default value for initial density of species
        call inp_read_array(u, 'initialspeciesdensity', inp%IniSpecDen)
        call inp_read_array(u, 'initialspeciesnumber',     inp%IniSpecNum)
        inp%SpecDenFix=0
        call inp_read_array(u, 'fixspeciesdensity',   inp%SpecDenFix)

        ! Diffusion events
        call inp_read_array(u, 'diffusionbarrier',      inp%DiffBrr, stat)
        if ( stat > 0 ) call print_err("Error in reading DiffusionBarrier")
    
        ! Att&Det event
        inp%attflag=0
        if (inp%FlagCore.eqv..True. ) then
            allocate(inp%AttBrr(inp%SpecNum,3),inp%DetBrr(inp%SpecNum,3))
            inp%AttBrr=-1
            inp%DetBrr=-1
            call inp_read_array(u, 'attachbarrier',     inp%AttBrr(:,1), stat)
            if ( stat > 0 ) then
                inp%attflag=0 !1 for all, 2 for only attachbarrier1 provided, 0 for nothing
            else 
                call inp_read_array(u, 'attachbarrier2',    inp%AttBrr(:,2), stat)
                if ( stat > 0 ) then
                    inp%AttBrr(:,2)=inp%AttBrr(:,1)
                    inp%AttBrr(:,3)=inp%AttBrr(:,1)
                    inp%attflag=2
                else
                    inp%attflag=1
                    call inp_read_array(u, 'attachbarrier3',    inp%AttBrr(:,3), stat)
                    if ( stat > 0 )  call print_err("attachbarrier3 should be provided")
                end if 
            end if
            if(inp%attflag/=0) then
                call inp_read_array(u, 'detachbarrier',     inp%DetBrr(:,1), stat)
                if(inp%attflag==1) then
                    call inp_read_array(u, 'detachbarrier2',    inp%DetBrr(:,2), stat)
                    if ( stat > 0 ) inp%DetBrr(:,2)=inp%DetBrr(:,1)
                    call inp_read_array(u, 'detachbarrier3',    inp%DetBrr(:,3), stat)
                    if ( stat > 0 ) inp%DetBrr(:,3)=inp%DetBrr(:,1)
                end if
            end if
        end if
    end if

  !-----------------------------Other events-------------------------------!
    ! Merge events
    call inp_read_value(u, 'mergeeventnumber',      inp%MergeNum, stat)
    if ( stat==0 .and. inp%MergeNum > 0) then
        allocate(inp%MergeR(inp%MergeNum,3),        inp%MergeBrr(inp%MergeNum))
        call inp_read_block(u, 'mergereaction',     inp%MergeR, inp%MergeBrr,stat)
        if ( stat > 0 ) call print_err("Error in reading MergeReaction")
        do i=1, inp%MergeNum
            if(inp%MergeBrr(i)<0) call print_err("Merge barrier cannot be negative!")
        end do
    end if

    ! Decomp events
    call inp_read_value(u, 'decompeventnumber',     inp%DecompNum, stat)
    if ( stat==0 .and. inp%DecompNum > 0) then
        allocate(inp%DecompR(inp%DecompNum,3),      inp%DecompBrr(inp%DecompNum))
        call inp_read_block(u, 'decompreaction',    inp%DecompR, inp%DecompBrr,  stat)
        if ( stat > 0 ) call print_err("Error in reading DecompReaction")
        do i=1,inp%DecompNum
            if(inp%DecompBrr(i)<0) call print_err("Decomposition barrier cannot be negative!")
        end do
    end if

    ! Desorp events
    call inp_read_value(u, 'desorpeventnumber',     inp%DesorpNum, stat)
    if ( stat==0 .and. inp%DesorpNum > 0) then
        allocate(inp%DesorpR(inp%DesorpNum),        inp%DesorpBrr(inp%DesorpNum))
        call inp_read_block(u, 'desorpreaction',    inp%DesorpR, inp%DesorpBrr,  stat)
        if ( stat > 0 ) call print_err("Error in reading DesorpReaction")
        do i=1,inp%DesorpNum
            if(inp%DesorpBrr(i)<0) call print_err("Desorption barrier cannot be negative!")
        end do
    end if

    close(u)

    !------------------------------------------------!

    ! for restart
    if(trim(io_lower(inp%Restart))=='true') then
        call teio_assert_file_exists("Status")
        u=io_unit()
        open(u,file="Status",status='old',action='read')

        call inp_read_value(u,"numberofspeciestype",inp%SpecNumR)
        ! SpecNumR can be smaller than SpecNum
        if(inp%SpecNumR>0) then
            allocate(inp%NumSpecies(inp%SpecNumR))
            call inp_read_array(u,"numberofeachspecies",inp%NumSpecies)
        end if

        call inp_read_value(u,"numberofmspeciestype",inp%MSpecNumR)
        ! MSpecNumR can be smaller than MSpecNum
        if(inp%MSpecNumR>0) then
            allocate(inp%NumMSpecies(inp%MSpecNumR))
            call inp_read_array(u,"numberofeachmspecies",inp%NumMSpecies)
        end if

        call inp_read_value(u,"gridoccnum",inp%GridOccNum)
        if(inp%GridOccNum>0) then
            allocate(inp%GridOcc(inp%GridOccNum,3))
            call inp_read_block(u,"GridOcc",inp%GridOcc)
        end if

        ! read Edge status
        call inp_read_array(u,"edgesitenum",inp%EdgeSiteNum,stat)
        if(stat==0) then
            if(inp%EdgeSiteNum(1)>0) then
                allocate(inp%Edge1(inp%EdgeSiteNum(1)))
                call inp_read_array(u,"edge1",inp%Edge1,stat)
                if(stat>0) inp%Edge1=0
            end if
            if(inp%EdgeSiteNum(2)>0) then
                allocate(inp%Edge2(inp%EdgeSiteNum(2)))
                call inp_read_array(u,"edge2",inp%Edge2,stat)
                if(stat>0) inp%Edge2=0
            end if
            if(inp%EdgeSiteNum(3)>0) then
                allocate(inp%Edge3(inp%EdgeSiteNum(3)))
                call inp_read_array(u,"edge3",inp%Edge3,stat)
                if(stat>0) inp%Edge3=0
            end if
            if(inp%EdgeSiteNum(4)>0) then
                allocate(inp%Edge4(inp%EdgeSiteNum(4)))
                call inp_read_array(u,"edge4",inp%Edge4,stat)
                if(stat>0) inp%Edge4=0
            end if
        end if

        close(u)
        ! check Status data
        if(sum(inp%NumSpecies)/=inp%GridOccNum) &
            call print_err("Error in analysing Status file content!")
    end if

    inp%init = .true.

  end function read_InpKMC

  !--------------------------------------------------------------------!
  !------------------------------- TE ---------------------------------!
  !--------------------------------------------------------------------!

  function read_InpTE(file) result(inp)

    implicit none

    character(len=*), intent(in) :: file
    type(InputData)              :: inp

    integer                :: iline, stat
    integer                :: u

    call teio_assert_file_exists(file)
    inp%file = trim(file)
    call inp_defaults(inp)

    u = io_unit()
    open(u, file=trim(adjustl(file)), status='old', action='read')

    close(u)

    inp%init = .true.

  end function read_InpTE

  !--------------------------------------------------------------------!
  !                         input data summary                         !
  !--------------------------------------------------------------------!

  subroutine inp_print_info(inp)
    implicit none

    type(InputData) :: inp
    if (.not. inp%init) return

    if (trim(io_lower(inp%mode))== 'kmc') then
       call inp_print_info_KMC(inp)
    else if (trim(io_lower(inp%mode))== 'te') then
       call inp_print_info_TE(inp)
    end if
  end subroutine inp_print_info

  !--------------------------------------------------------------------!

  subroutine inp_print_info_KMC(inp)
    implicit none

    type(InputData) :: inp
    integer         :: i

    if(trim(io_lower(inp%Restart))=='true') then
        write(*,'(1x,70("="))')
        write(*,*) 'Restart simulation'
    end if
    write(*,'(1x,70("="))')
    write(*,*)
    write(*,*) 'Kinetic Monte-Carlo Simulation Parameters'
    write(*,*)
    write(*,*) 'Temperature : ', trim(io_adjustl(inp%T)), ' K'
    write(*,"(A,E11.4,A)") ' Vibration frequency : ', inp%Vib, ' s^-1'

    if(inp%MeshType==2) then
        write(*,*) 'Mesh type   : Square'
    else if(inp%MeshType==1) then
        write(*,*) 'Mesh type   : Hexagonal'
    end if
    write(*,*) 'Mesh size   : ', trim(io_adjustl(inp%MeshSize(1))),' x ',trim(io_adjustl(inp%MeshSize(2)))
    if(inp%FlagCore)  write(*,*) 'Core size   : ', trim(io_adjustl(inp%CoreSize(1))),' x ',trim(io_adjustl(inp%CoreSize(2)))
    write(*,*)

    write(*,*) 'Total step                               : ', trim(io_adjustl(inp%TotStep))
    write(*,*) 'Step interval for species density print  : ', trim(io_adjustl(inp%PrtSpecDenStep))
    write(*,*) 'Step interval for current status save    : ', trim(io_adjustl(inp%SaveStatusStep))
    write(*,*)


    if(inp%GasNum>0) then     ! Gas
        write(*,*) 'Gas type number : ', trim(io_adjustl(inp%GasNum))
        write(*,*) '---------------------------------------------------'
        write(*,*) '|   Gas   |  Partial  | Adsorption |  Adsorption  | ' 
        write(*,*) '|   Name  |  pressure |    rate    |   product    | '
        write(*,*) '---------------------------------------------------'
        do i=1,inp%GasNum
          write(*,"(' |',A7$)")  trim(adjustl(inp%GasName(i)))
          write(*,"('  |',E10.3$)")  inp%GasPrs(i)
          write(*,"(' |',E11.4,' '$)")  inp%GasAdsorbRate(i)
          if(trim(adjustl(inp%GasProd(i,3)))=="") then
            write(*,"('|  ',A7,'     |')") trim(adjustl(inp%GasProd(i,2)))
          else 
            write(*,"('|  ',A4,A2,A4,'  |')") trim(adjustl(inp%GasProd(i,2)))," +",trim(adjustl(inp%GasProd(i,3)))
          end if
        end do
        write(*,*) '---------------------------------------------------'
        write(*,*)
    end if

    if(inp%MSpecNum+inp%SpecNum>0) then
        write(*,*) 'Species type number :  ', trim(io_adjustl(inp%MSpecNum+inp%SpecNum))
        if (inp%FlagCore) then
            if(inp%attflag==1.or.inp%mattflag==1) then
                write(*,'(A)') ' -------------------------------------------------------------------------------------'
                write(*,'(A)') ' | Species | MFA | Initial |Diffusion|  Att  |  Att2 |  Att3 |  Det  |  Det2 |  Det3 |' 
                write(*,'(A)') ' |  Name   |     | Density | Barrier |Barrier|Barrier|Barrier|Barrier|Barrier|Barrier|'
                write(*,'(A)') ' -------------------------------------------------------------------------------------'
            else 
                write(*,'(A)') ' -----------------------------------------------------'
                write(*,'(A)') ' | Species | MFA | Initial |Diffusion|Attach |Dettach|' 
                write(*,'(A)') ' |  Name   |     | Density | Barrier |Barrier|Barrier|'
                write(*,'(A)') ' -----------------------------------------------------'
            end if
        else 
            write(*,'(A)') ' -------------------------------------'
            write(*,'(A)') ' | Species | MFA | Initial |Diffusion|' 
            write(*,'(A)') ' |  Name   |     | Density | Barrier |'
            write(*,'(A)') ' -------------------------------------'
        end if
    end if

    if(inp%MSpecNum>0) then
        do i = 1, inp%MSpecNum
            write(*,"(' |',A7$)")  trim(adjustl(inp%MSpecName(i)))
            write(*,"('  | Yes '$)")  !MFA or not
            if(inp%MIniSpecDen(i)>0) then ! write initial density
                write(*,"('| ',F7.5$)")  inp%MIniSpecDen(i)
            else if(inp%MIniSpecNum(i)>0) then
                write(*,"('| ',F7.5$)") &  
                    real(inp%MIniSpecNum(i))/(inp%MeshSize(1)*inp%MeshSize(2)-inp%CoreSize(1)*inp%CoreSize(2))
            else 
                write(*,"('| ',F7.5$)")  0.0
            end if
            write(*,"(' |   ---  '$)")  ! diffusion
            if(inp%FlagCore.and.inp%mattflag/=0) then
                write(*,"(' |',F6.3$)")      inp%MAttBrr(i,1)
                if(inp%mattflag==1) then
                    write(*,"(' |',F6.3$)")  inp%MAttBrr(i,2)
                    write(*,"(' |',F6.3$)")  inp%MAttBrr(i,3)
                end if
                write(*,"(' |',F6.3$)")      inp%MDetBrr(i,1)
                if(inp%mattflag==1) then
                    write(*,"(' |',F6.3$)")  inp%MDetBrr(i,2)
                    write(*,"(' |',F6.3$)")  inp%MDetBrr(i,3)
                end if
                write(*,"(' |')")
            else
                write(*,"(' |')")
            end if
        end do
        if (inp%FlagCore) then
            if(inp%attflag==1) then
                write(*,'(A)') ' ------------------------------------------------------------------------------------- '
            else
                write(*,'(A)') ' ----------------------------------------------------- '
            end if
        else
            write(*,'(A)') ' ------------------------------------- '
        end if
    end if

    if(inp%SpecNum>0) then    ! Species
        do i = 1, inp%SpecNum
            write(*,"(' |',A7$)")  trim(adjustl(inp%SpecName(i)))
            write(*,"('  | No'$)") 
            if(inp%IniSpecDen(i)>0) then 
                write(*,"(' | ',F7.5$)")  inp%IniSpecDen(i)
            else if(inp%IniSpecNum(i)>0) then
                write(*,"('  | ',F7.5$)") &
                    real(inp%IniSpecNum(i))/(inp%MeshSize(1)*inp%MeshSize(2)-inp%CoreSize(1)*inp%CoreSize(2))
            else 
                write(*,"('  | ',F7.5$)")  0.0
            end if
            !write(*,"(' |    ---   '$)") 
            if(inp%DiffBrr(i)>0) then 
                write(*,"(' |',F7.3,' '$)")  inp%DiffBrr(i)
            else 
                write(*,"(' |   ---  '$)") 
            end if
            if(inp%FlagCore.and.inp%attflag/=0) then
                write(*,"(' |',F6.3$)")      inp%AttBrr(i,1)
                if(inp%attflag==1) then
                    write(*,"(' |',F6.3$)")  inp%AttBrr(i,2)
                    write(*,"(' |',F6.3$)")  inp%AttBrr(i,3)
                end if
                write(*,"(' |',F6.3$)")      inp%DetBrr(i,1)
                if(inp%attflag==1) then
                    write(*,"(' |',F6.3$)")  inp%DetBrr(i,2)
                    write(*,"(' |',F6.3$)")  inp%DetBrr(i,3)
                end if
                write(*,"(' |')")
            else
                write(*,"(' |')")
            end if
        end do
    end if

    if(inp%MSpecNum+inp%SpecNum>0) then
        if (inp%FlagCore) then
            if(inp%attflag==1) then
                write(*,'(A)') ' ------------------------------------------------------------------------------------- '
            else
                write(*,'(A)') ' ----------------------------------------------------- '
            end if
        else
            write(*,'(A)') ' ------------------------------------- '
        end if
        write(*,'(A)') 
    end if

    if(inp%MSpecNum>0) then
       do i=1, inp%MSpecNum
          if(inp%MSpecDenFix(i)/=0) then
            write(*,*) "Density of MSpecies: "//trim(adjustl(inp%MSpecName(i)))//" is fixed."
          end if
       end do
    end if
    write(*,*) 

    if(inp%SpecNum>0) then  
       do i=1, inp%SpecNum
          if(inp%SpecDenFix(i)/=0) then
            write(*,*) "Density of Species: "//trim(adjustl(inp%SpecName(i)))//" is fixed."
          end if
       end do
    end if
    write(*,*) 

10  Format(A6,A,A6,A,A6,A,F7.3)       
    if(inp%MergeNum>0) then  ! Merge event info.
        write(*,*) 'Merge event number : ', trim(io_adjustl(inp%MergeNum))
        !write(*,*) '---------------------------------------------------'
        do i=1,inp%MergeNum 
            write(*,10) trim(adjustl(inp%MergeR(i,1))),'   +',trim(adjustl(inp%MergeR(i,2))),&
            '   ->',trim(adjustl(inp%MergeR(i,3))),'   with barrier : ',inp%MergeBrr(i)
        end do
        !write(*,*) '---------------------------------------------------'
        write(*,*)
    end if

    if(inp%DecompNum>0) then ! Decomposition event info.
        write(*,*) 'Decomposition event number : ', trim(io_adjustl(inp%DecompNum))
        !write(*,*) '---------------------------------------------------'
        do i=1,inp%DecompNum 
            write(*,10) trim(adjustl(inp%DecompR(i,1))),'   ->',trim(adjustl(inp%DecompR(i,2))),&
            '   +',trim(adjustl(inp%DecompR(i,3))),'   with barrier : ',inp%DecompBrr(i)
        end do
        !write(*,*) '---------------------------------------------------'
        write(*,*)
    end if

    if(inp%DesorpNum>0) then ! Desorption event info.
        write(*,*) 'Desorption event number : ', trim(io_adjustl(inp%DesorpNum))
        !write(*,*) '---------------------------------------------------'
        do i=1,inp%DesorpNum 
            write(*,'(A6,A,F7.3)') trim(adjustl(inp%DesorpR(i))),' desorption with barrier : ',inp%DesorpBrr(i)
        end do
        !write(*,*) '---------------------------------------------------'
        write(*,*)
    end if

    write(*,'(1x,70("="))')
    write(*,*)

  end subroutine inp_print_info_KMC

  subroutine inp_print_info_TE(inp)
    implicit none
    type(InputData) :: inp
    integer         :: i
    write(*,'(1x,70("="))')
  
  end subroutine inp_print_info_TE
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !     simply read value(s) for specific keyword (implementation)     !
  !--------------------------------------------------------------------!

  subroutine inp_read_value_c1(unit, keyword, dest, stat)
    implicit none
    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    character(len=*),  intent(inout) :: dest
    integer, optional, intent(out)   :: stat
    integer                          :: iline, ipos
    character(len=LINELEN)           :: line
    character(len=10)                :: str
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, dest)
    end if
  end subroutine inp_read_value_c1

  !--------------------------------------------------------------------!

  subroutine inp_read_value_i1(unit, keyword, dest, stat)
    implicit none
    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    integer,           intent(inout) :: dest
    integer, optional, intent(out)   :: stat
    integer                          :: iline, ipos
    character(len=LINELEN)           :: line
    character(len=10)                :: str
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, dest)
    end if
  end subroutine inp_read_value_i1

  !--------------------------------------------------------------------!

  subroutine inp_read_value_i2(unit, keyword, dest, stat)
    implicit none
    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    integer(kind=8),   intent(inout) :: dest
    integer, optional, intent(out)   :: stat
    integer                          :: iline, ipos
    character(len=LINELEN)           :: line
    character(len=10)                :: str
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, dest)
    end if
  end subroutine inp_read_value_i2

  !--------------------------------------------------------------------!
  
  subroutine inp_read_value_d1(unit, keyword, dest, stat)
    implicit none
    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    double precision,  intent(inout) :: dest
    integer, optional, intent(out)   :: stat
    integer                :: iline, ipos
    character(len=LINELEN) :: line
    character(len=10)      :: str
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, dest)
    end if
  end subroutine inp_read_value_d1

  !--------------------------------------------------------------------!
  ! Flags can be specified in three different ways in the input file   !
  !                                                                    !
  !   (1) KEYWORD                                                      !
  !   (2) KEYWORD .true.                                               !
  !   (3) KEYWORD .false.                                              !
  !                                                                    !
  ! Both, (1) and (2) will result in the return value '.true.', while  !
  ! (3) will result in '.false.'.  If KEYWORD is not found, the input  !
  ! value of flag will be returned.                                    !
  !--------------------------------------------------------------------!

  subroutine inp_read_value_l1(unit, keyword, flag, stat)
    implicit none
    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    logical,           intent(inout) :: flag
    integer, optional, intent(out)   :: stat
    integer                :: iline, ipos
    character(len=LINELEN) :: line
    character(len=10)      :: str
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, str)
       if (ipos > 0) then
          read(str,*) flag
       else
          flag = .true.
       end if
    end if
  end subroutine inp_read_value_l1

  !--------------------------------------------------------------------!
  !                  procedures for specific keywords                  !
  !--------------------------------------------------------------------!

  subroutine inp_read_array_c(unit, keyword, dest, stat)
    implicit none
    integer,                              intent(in)    :: unit
    character(len=*),                     intent(in)    :: keyword
    character(len=SPCSLEN), dimension(:), intent(inout) :: dest
    integer, optional,                    intent(out)   :: stat
    integer :: i, iline, ipos, asize
    character(len=LINELEN) :: line
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline)
    asize=size(dest)
    if (iline > 0) then
       call teio_readline(unit, iline, line)
       ipos=1
       do i=1,asize
           call io_readnext(line, ipos, dest(i))
       end do
    else
       if (present(stat)) stat = S_NOT
    end if
  end subroutine inp_read_array_c

  !--------------------------------------------------------------------!

  subroutine inp_read_array_i(unit, keyword, dest, stat)
    implicit none
    integer,                intent(in)    :: unit
    character(len=*),       intent(in)    :: keyword
    integer, dimension(:),  intent(inout) :: dest
    integer, optional,      intent(out)   :: stat
    integer :: i, iline, ipos, asize
    character(len=LINELEN) :: line
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline)
    asize=size(dest)
    if (iline > 0) then
       call teio_readline(unit, iline, line)
       ipos=1
       do i=1,asize
           call io_readnext(line, ipos, dest(i))
       end do
    else
       if (present(stat)) stat = S_NOT
    end if
  end subroutine inp_read_array_i

  !--------------------------------------------------------------------!

  subroutine inp_read_array_d(unit, keyword, dest, stat)
    implicit none
    integer,                         intent(in)    :: unit
    character(len=*),                intent(in)    :: keyword
    double precision,  dimension(:), intent(inout) :: dest
    integer, optional,               intent(out)   :: stat
    integer :: i, iline, ipos, asize
    character(len=LINELEN) :: line
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline)
    asize=size(dest)
    if (iline > 0) then
       call teio_readline(unit, iline, line)
       ipos=1
       do i=1,asize
           call io_readnext(line, ipos, dest(i))
       end do
    else
       if (present(stat)) stat = S_NOT
    end if
  end subroutine inp_read_array_d

  !---------------------------------------------------------------------!

  subroutine inp_read_block_im(unit, keyword, dest1,dest2,stat) 
    ! for integer matrix dest1
    implicit none
    integer,                                intent(in)    :: unit
    character(len=*),                       intent(in)    :: keyword
    integer,                dimension(:,:), intent(inout) :: dest1
    double precision,optional,dimension(:), intent(inout) :: dest2
    integer, optional,                      intent(out)   :: stat
    integer :: i, j, asize, bsize, csize
    integer :: iline, ipos
    character(len=LINELEN) :: line
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline)
    asize=size(dest1,dim=1)
    bsize=size(dest1,dim=2)
    if(present(dest2)) then
      csize=size(dest2)
      if (asize/=csize) then
         stat = S_ERROR
         return
      end if
    end if
    if (iline > 0) then
       do i = 1, asize
           call teio_readline(unit, iline, line)
           ipos=1
           do j=1,bsize
               call io_readnext(line, ipos, dest1(i,j))
           end do
           if(present(dest2)) call io_readnext(line, ipos, dest2(i))
       end do
    else
        if (present(stat)) stat = S_NOT
    end if
  end subroutine inp_read_block_im

  !--------------------------------------------------------------------!

  subroutine inp_read_block_cm(unit, keyword, dest1,dest2,stat) 
    ! for chars matrix dest1
    implicit none
    integer,                                intent(in)    :: unit
    character(len=*),                       intent(in)    :: keyword
    character(len=SPCSLEN), dimension(:,:), intent(inout) :: dest1
    double precision,optional,dimension(:), intent(inout) :: dest2
    integer, optional,                      intent(out)   :: stat
    integer :: i, j, asize, bsize, csize
    integer :: iline, ipos
    character(len=LINELEN) :: line
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline)
    asize=size(dest1,dim=1)
    bsize=size(dest1,dim=2)
    if(present(dest2)) then
      csize=size(dest2)
      if (asize/=csize) then
         stat = S_ERROR
         return
      end if
    end if
    if (iline > 0) then
       do i = 1, asize
           call teio_readline(unit, iline, line)
           ipos=1
           do j=1,bsize
               call io_readnext(line, ipos, dest1(i,j))
           end do
           if(present(dest2)) call io_readnext(line, ipos, dest2(i))
       end do
    else
        if (present(stat)) stat = S_NOT
    end if
  end subroutine inp_read_block_cm
  !--------------------------------------------------------------------!

  subroutine inp_read_block_ca(unit, keyword, dest1,dest2,stat) 
    ! for chars array dest1
    implicit none
    integer,                                intent(in)    :: unit
    character(len=*),                       intent(in)    :: keyword
    character(len=SPCSLEN), dimension(:),   intent(inout) :: dest1
    double precision,optional,dimension(:), intent(inout) :: dest2
    integer, optional,                      intent(out)   :: stat
    integer :: i, asize, csize
    integer :: iline, ipos
    character(len=LINELEN) :: line
    if (present(stat)) stat = S_OK
    call inp_find_keyword(unit, keyword, iline)
    asize=size(dest1)
    if(present(dest2)) then
      csize=size(dest2)
      if (asize/=csize) then
         stat = S_ERROR
         return
      end if
    end if
    if (iline > 0) then
       do i = 1, asize
           call teio_readline(unit, iline, line)
           ipos=1
           call io_readnext(line, ipos, dest1(i))
           if(present(dest2)) call io_readnext(line, ipos, dest2(i))
       end do
    else
        if (present(stat)) stat = S_NOT
    end if
  end subroutine inp_read_block_ca

  !--------------------------------------------------------------------!
  !                         read until keyword                         !
  !--------------------------------------------------------------------!

  subroutine inp_find_keyword(unit, keyword, iline, line)
    implicit none
    integer,                          intent(in)  :: unit
    character(len=*),                 intent(in)  :: keyword
    integer,                          intent(out) :: iline
    character(len=LINELEN), optional, intent(out) :: line
    integer                :: stat
    character(len=LINELEN) :: kwd, ln
    rewind(unit)
    kwd = ''
    iline = 0
    skip : do
       call teio_readline(unit, iline, ln, stat)
       if (stat /= 0) then
          ! end of file and keyword not found
          iline = 0
          exit skip
       end if
       iline = iline + 1
       read(ln, *) kwd
       if (trim(io_lower(kwd))==trim(io_lower(keyword))) then
          ! keyword found
          if (present(line)) line = ln
          exit skip
       end if
    end do skip
  end subroutine inp_find_keyword

  !--------------------------------------------------------------------!

  subroutine print_err(err_info,eflag)
    character(len=*),          intent(in)  ::  err_info
    integer, optional,         intent(in)  ::  eflag
    if(.not.present(eflag) .or. eflag==1) &
        write(*,*) "------------------Error!!!------------------"
    write(*,*) err_info
    if(.not.present(eflag) .or. eflag==3) then
        write(*,*) "------------------Error!!!------------------"
        stop
    end if
  end subroutine print_err

  !--------------------------------------------------------------------!
  
  function read_long_int(LongF) result(LongInt)
    double precision,    intent(in)  :: LongF
    integer(kind=8)                  :: LongInt
    character(19)                    :: chars
    integer                          :: i
    write(chars,'(f19.0)') LongF
    i=scan(chars,'.')-1
    read(chars(1:i),*) LongInt
  end function read_long_int
  
end module input
