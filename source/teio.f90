!-----------------------------------------------------------------------
!                           teio.f90 
!-----------------------------------------------------------------------

module teio

  use io, only: io_lower,    &
                io_readval,  &
                io_adjustl,  &
                io_center
  implicit none

  public :: teio_readline,             &
            teio_timestamp,            &
            teio_assert_file_exists

  private :: teio_readline_c1, teio_readline_cn, &
             teio_readline_i1, teio_readline_in, &
             teio_readline_d1, teio_readline_dn, &
             teio_readline_l1, teio_readline_ln

  !---------------------------- constants -----------------------------!
  ! LINELEN    max. length of a line read from an input file           !
  ! PATHLEN    max. number of characters available for a file path     !
  ! TYPELEN    max. number of characters for atmoic species names      !
  ! STDIN      file unit of standard in                                !
  ! STDOUT     file unit of standard out                               !
  ! STDERR     file unit of standard error                             !
  !--------------------------------------------------------------------!

  integer, parameter, public :: LINELEN = 1024
  integer, parameter, public :: PATHLEN = 1024
  integer, parameter, public :: SPCSLEN = 10
  integer, parameter, public :: TYPELEN = 2
  integer, parameter, public :: STDIN   = 5
  integer, parameter, public :: STDOUT  = 6
  integer, parameter, public :: STDERR  = 0

  !--------------------------------------------------------------------!
  !   teio_readline() - read next line with contents from input file   !
  !                                                                    !
  ! A line is skipped if                                               !
  !                                                                    !
  ! - it only contains blanks                                          !
  ! - the first non-blank character is `!', `#', or `%'                !
  !                                                                    !
  ! usage: call teio_readline(unit, iline, line[, n][, stat])          !
  !        implementations available for 'line' as character, integer, !
  !        double precision, and logical (also array of length 'n')    !
  !--------------------------------------------------------------------!

  interface teio_readline
     module procedure teio_readline_c1, teio_readline_cn, &
                      teio_readline_i1, teio_readline_in, &
                      teio_readline_d1, teio_readline_dn, &
                      teio_readline_l1, teio_readline_ln
  end interface teio_readline

contains

  !--------------------------------------------------------------------!
  !              return a formatted date and time string               !
  !--------------------------------------------------------------------!

  function teio_timestamp() result(date)

    implicit none

    character(len=30)     :: date
    integer, dimension(8) :: v

    call date_and_time(values=v)

    write(date, '(I4.4,"-",I2.2,"-",I2.2,2x,I2.2,":",I2.2,":",I2.2)') &
         v(1:3), v(5:7)

  end function teio_timestamp

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine teio_assert_file_exists(file)
    implicit none
    character(len=*), intent(in) :: file
    logical :: fexists
    inquire(file=trim(adjustl(file)), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: file not found: ", trim(adjustl(file))
       stop
    end if
  end subroutine teio_assert_file_exists

  !--------------------------------------------------------------------!
  !     Implementation of teio_readline() for different data types     !
  !--------------------------------------------------------------------!

  subroutine teio_readline_c1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    character(len=*),  intent(out)   :: line
    integer, optional, intent(out)   :: stat

    integer :: stat2

    stat2 = 0
    do
       read(u_in, '(A)', iostat=stat2) line
       if (stat2 == 0) then
          iline = iline + 1
          line  = trim(adjustl(line))
          if (line(1:1) == '!')    cycle
          if (line(1:1) == '#')    cycle
          if (line(1:1) == '%')    cycle
          if (len_trim(line) == 0) cycle
       end if
       exit
    end do
    if (present(stat)) stat = stat2

  end subroutine teio_readline_c1

  !--------------------------------------------------------------------!

  subroutine teio_readline_cn(u_in, iline, line, n, stat)

    implicit none

    integer,                        intent(in)    :: u_in
    integer,                        intent(inout) :: iline
    integer,                        intent(in)    :: n
    character(len=*), dimension(n), intent(out)   :: line
    integer, optional,              intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call teio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine teio_readline_cn

  !--------------------------------------------------------------------!

  subroutine teio_readline_i1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    integer,           intent(out)   :: line
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call teio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line
    end if
    if (present(stat)) stat = stat2

  end subroutine teio_readline_i1

  !--------------------------------------------------------------------!

  subroutine teio_readline_in(u_in, iline, line, n, stat)

    implicit none

    integer,               intent(in)    :: u_in
    integer,               intent(inout) :: iline
    integer,               intent(in)    :: n
    integer, dimension(n), intent(out)   :: line
    integer, optional,     intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call teio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine teio_readline_in

  !--------------------------------------------------------------------!

  subroutine teio_readline_d1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    double precision,  intent(out)   :: line
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call teio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line
    end if
    if (present(stat)) stat = stat2

  end subroutine teio_readline_d1

  !--------------------------------------------------------------------!

  subroutine teio_readline_dn(u_in, iline, line, n, stat)

    implicit none

    integer,                        intent(in)    :: u_in
    integer,                        intent(inout) :: iline
    integer,                        intent(in)    :: n
    double precision, dimension(n), intent(out)   :: line
    integer, optional,              intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call teio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine teio_readline_dn

  !--------------------------------------------------------------------!

  subroutine teio_readline_l1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    logical,           intent(out)   :: line
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call teio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line
    end if
    if (present(stat)) stat = stat2

  end subroutine teio_readline_l1

  !--------------------------------------------------------------------!

  subroutine teio_readline_ln(u_in, iline, line, n, stat)

    implicit none

    integer,                        intent(in)    :: u_in
    integer,                        intent(inout) :: iline
    integer,                        intent(in)    :: n
    logical, dimension(n),          intent(out)   :: line
    integer, optional,              intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call teio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine teio_readline_ln

end module teio
