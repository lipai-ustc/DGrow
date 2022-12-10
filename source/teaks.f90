!-----------------------------------------------------------------------!
!  This is the main program of the TEAKS package.                       !
!-----------------------------------------------------------------------!
!  Copyright (C) 2020-2026 Pai Li lipai@mail.ustc.edu.cn                !
!  2020-01-27 Pai Li                                                    !
!-----------------------------------------------------------------------!

program teaks

  use teio,     only: teio_readline,        &
                      teio_timestamp,       &
                      PATHLEN, LINELEN, SPCSLEN

  use input,    only: InputData,      &
                      read_InpKMC,    &
                      inp_print_info

  use io,       only: io_adjustl,     &
                      io_center,      &
                      io_unit

  use mesh,     only: init_mesh,      &
                      save_grid,      &
                      SiteNum

  use core,     only: save_core,      &
                      EdgeProb,       &
                      init_core,      &
                      Edges

  use species,  only: init_species,   &
                      SpecProb,       &
                      MSpecProb,      &
                      MMMergeProb,    &
                      species_step,   &
                      mspecies_step,  &
                      mmmerge_step,   &
                      detach,         &
                      save_species,   &
                      print_species_info,&
                      print_species_density,&
                      print_species_name,&
                      species_prob_update,&
                      fix_den

  use gas,      only: init_gas,       &
                      GasProb,        &
                      adsorption,     &
                      print_gas_info

  use statistics, only: init_statistics,   &
                        print_event_count

  implicit none

  !--------------------------------------------------------------------

  type(InputData)                :: inp
  character(len=PATHLEN)         :: inFile

  integer(kind=8)                :: istep
  double precision               :: ProbSum
  double precision               :: ProbRand
  double precision               :: time
  double precision               :: randn
  integer                        :: u_num,u_den,u_edge

  character(len=30)  :: time_start

  !-------------------------- initialization --------------------------!
  time_start=teio_timestamp()
  call random_init()
  call initialize(inFile)

  inp = read_InpKMC(inFile)
  if ( .not. inp%init ) stop
  call inp_print_info(inp)

  call init_statistics(inp)
  ! model initialization. first mesh,then species, finally gas
  call init_mesh(inp)
  call init_species(inp)
  call init_gas(inp)
  if(inp%FlagCore) call init_core(inp)
  call print_event_info()

  write(*,*) "---------------------------------------------------------------------------"
  write(*,*) "    istep  ProbSum    GasProb    SpecProb   MSpecProb MMMergeProb EdgeProb"
  write(*,*) "---------------------------------------------------------------------------"

  u_num=88
  u_den=99
  u_edge=98
  open(u_num,file='SpecNum',status='replace',action='write')
  open(u_den,file='Density',status='replace',action='write')
  open(u_edge,file='EdgeInfo',status='replace',action='write')
  write(u_num,*) "Site Number: ", trim(io_adjustl(SiteNum))
  write(u_num,'(A$)') " Time (s)  "
  write(u_den,'(A$)') " Time (s)  "
  write(u_edge,'(A,4(I4," "))') " EdgeSiteNum: ",Edges(1)%enum,Edges(2)%enum,Edges(3)%enum,Edges(4)%enum
  write(u_edge,'(A)') " EdgeOccNum:"
  call print_species_name(u_num)
  call print_species_name(u_den)
  !---------------------------- simulation ----------------------------!
  time=0 !unit in s

  write(u_num,'(E10.3,x$)') time
  write(u_den,'(E10.3,x$)') time
  call print_species_density(u_num,1)
  call print_species_density(u_den)
  ProbSum=GasProb+SpecProb+MSpecProb+MMMergeProb+EdgeProb
  write(*,'(I10,6(x,E10.3))') 0,ProbSum,GasProb,SpecProb,MSpecProb,MMMergeProb,EdgeProb
      
  do istep=1, inp%TotStep
      !write(*,*) "Step: ",istep
      ProbSum=GasProb+SpecProb+MSpecProb+MMMergeProb+EdgeProb
      call random_number(randn)
      ProbRand=ProbSum*randn
      if(GasProb>ProbRand) then
          !write(*,*) "teaks adsorption"
          call adsorption()
      else if(GasProb+SpecProb>ProbRand) then
          !write(*,*) "teaks species_step"
          call species_step()
      else if(GasProb+SpecProb+MSpecProb>ProbRand) then
          !write(*,*) "teaks mspecies_step"
          call mspecies_step()
      else if(ProbSum-EdgeProb>ProbRand) then
          !write(*,*) "teaks mmmerge_step"
          call mmmerge_step()
      else
          call detach()
      end if

      !do 
          call random_number(randn)
      !    if(randn/=0) exit
      !end do
      time=time+log(1/randn)/ProbSum

      if(mod(istep,1000)==0) call fix_den()

      if(mod(istep,inp%PrtSpecDenStep)==0)  then
          write(u_num,'(E10.3,x$)') time
          call print_species_density(u_num,1)
          write(u_den,'(E10.3,x$)') time
          call print_species_density(u_den)
      end if

      if(mod(istep,inp%SaveStatusStep)==0) then
          call save_status()
          write(*,'(I10,6(x,E10.3))') istep,ProbSum,GasProb,SpecProb,MSpecProb,MMMergeProb,EdgeProb
          write(u_edge,'(4(I4," "))') Edges(1)%snum,Edges(2)%snum,Edges(3)%snum,Edges(4)%snum
      end if
  end do

  call species_prob_update()
  write(*,'(I10,6(x,E10.3))') istep,ProbSum,GasProb,SpecProb,MSpecProb,MMMergeProb,EdgeProb
  write(*,*) "---------------------------------------------------------------------------"

  close(u_num)
  close(u_den)
  close(u_edge)

  call finalize(.True.)

contains !=============================================================!

  subroutine initialize(inFile)
    implicit none
    character(len=*), intent(out) :: inFile
    integer :: nargs
    logical :: fexists

    write(*,*) repeat('=',70)
    write(*,*) io_center("Thermodynamic equilibrium analysis and KMC simulation (TEAKS)", 70)
    write(*,*) repeat('=',70)
    write(*,*)
    write(*,*) "Copyright (C) 2020-2026 Pai Li @ USTC"

    nargs = command_argument_count()
    if (nargs < 1) then
       write(0,*) "Error: No input file provided."
       call print_usage()
       call finalize(.False.)
       stop
    end if

    call get_command_argument(1, value=inFile)
    inquire(file=trim(inFile), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: File not found: ", trim(inFile)
       call print_usage()
       call finalize(.False.)
       stop
    end if

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine save_status()
    implicit none
    integer       :: u
    u=io_unit()
    open(u,file='Status',status='replace',action='write')
    call save_species(u)
    call save_grid(u)
    call save_core(u)
    close(u)
  end subroutine save_status

  !--------------------------------------------------------------------!

  subroutine print_event_info()
    implicit none
    integer       :: u
    u=io_unit()
    open(u,file='EventInfo',status='replace',action='write')
    call print_gas_info(u)
    call print_species_info(u)
    close(u)
  end subroutine print_event_info

  !--------------------------------------------------------------------!

  subroutine print_statistics()
    implicit none
    integer        :: u,i
    u=io_unit()
    open(u,file='Statistics',status='replace',action='write')
    write(u,'(A,x,E12.5,x,A)') "Time: ",time," s"
    write(u,'(A)') "Species Density: "
    call print_species_name(u)
    call print_species_density(u,1)
    write(u,*)
    write(u,'(A)') "==============Events Count============== "
    call print_event_count(u)
    if(inp%FlagCore)  then
        write(u,'(A$)')  "Edge front: "
        do i=1,4
            write(u,'(A," "$)')  trim(io_adjustl(Edges(i)%front))
        end do
        write(u,*)
    end if
    close(u)
  end subroutine

  !--------------------------------------------------------------------!

  subroutine finalize(fin)
    implicit none
    logical, intent(in) :: fin
    write(*,*)
    write(*,*) repeat('=',70)
    if(fin) then
        call save_status()
        call print_statistics()
        write(*,*) io_center("Program finished successfully",70)
    else
        write(*,*) io_center("Program finished unsuccessfully",70)
    end if
    write(*,*) io_center("Start at:  "//time_start,70)
    write(*,*) io_center("End   at:  "//teio_timestamp(),70)
    write(*,*) repeat('=',70)
    write(*,*)
  end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

    write(*,*)
    write(*,*) "teaks.x -- Thermodynamic equilibirum analysis and KMC simulation (TEAKS) "
    write(*,'(1x,70("-"))')
    write(*,*) 'Usage: teaks.x <input-file>'
    write(*,*)
    write(*,*) 'See the documentation or the source code for a description of the '
    write(*,*) 'input file format.'
    write(*,*)

  end subroutine print_usage

  !--------------------------------------------------------------------!
  subroutine random_init()
    implicit none
    integer :: i, n
    integer :: clock
    integer, dimension(:), allocatable :: iSeed
    i = 0 
    call random_seed(size=n)
    allocate(iSeed(n))
    call system_clock(count=clock)
    iSeed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=iSeed)
  end subroutine random_init

end program teaks
