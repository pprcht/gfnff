! ──────────────────────────────────────────────────────────────────────────────
! This file is part of gfnff.
!
! Copyright (C) 2023-2026 Philipp Pracht
!
! gfnff is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! gfnff is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with gfnff. If not, see <https://www.gnu.org/licenses/>.
! ──────────────────────────────────────────────────────────────────────────────
module gfnff_timer_mod
  !> Provides a simple wall-clock timer type for collecting labelled timings.
  use iso_fortran_env, only: real64, int64, output_unit
  implicit none
  private
  public :: gfnff_timer

  integer, parameter :: wp = real64
  integer, parameter :: label_len = 64

  type :: gfnff_timer
    !> Number of timing slots.
    integer :: n = 0
    !> System-clock tick at start for each slot (-1 = not started).
    integer(int64), allocatable :: t_start(:)
    !> Accumulated elapsed time in seconds for each slot (-1 = not finished).
    real(wp), allocatable :: t_elapsed(:)
    !> Label associated with each slot.
    character(len=label_len), allocatable :: labels(:)
    !> System clock ticks per second.
    integer(int64) :: count_rate = 0_int64
  contains
    procedure :: new     => timer_new
    procedure :: measure => timer_measure
    procedure :: write   => timer_write
  end type gfnff_timer

contains

  subroutine timer_new(self, n)
    !***********************************************
    !* Allocate space for n timing slots and       *
    !* initialise the timer object.                *
    !***********************************************
    class(gfnff_timer), intent(inout) :: self
    integer,            intent(in)    :: n
    !> n : number of independent timings to store

    integer(int64) :: dummy_count, dummy_max

    self%n = n

    if (allocated(self%t_start))   deallocate(self%t_start)
    if (allocated(self%t_elapsed)) deallocate(self%t_elapsed)
    if (allocated(self%labels))    deallocate(self%labels)

    allocate(self%t_start(n),   source=-1_int64)
    allocate(self%t_elapsed(n), source=-1.0_wp)
    allocate(self%labels(n))
    self%labels = ''

    call system_clock(count=dummy_count, count_rate=self%count_rate, &
         &            count_max=dummy_max)

  end subroutine timer_new


  subroutine timer_measure(self, i, label)
    !***********************************************
    !* Start or stop a timing slot.               *
    !*                                             *
    !* If label is present, the timer for slot i  *
    !* is started and the label recorded.         *
    !* If label is absent, the timer for slot i   *
    !* is stopped and the elapsed time stored.    *
    !*                                             *
    !* Input:                                      *
    !*   i     – slot index (1 .. n)              *
    !*   label – (optional) description string    *
    !***********************************************
    class(gfnff_timer), intent(inout)        :: self
    integer,            intent(in)           :: i
    character(len=*),   intent(in), optional :: label

    integer(int64) :: ticks

    if (i < 1 .or. i > self%n) return

    call system_clock(count=ticks)

    if (present(label)) then
      ! ---- start ----
      self%t_start(i) = ticks
      self%labels(i)  = label
    else
      ! ---- stop ----
      if (self%t_start(i) < 0_int64) return          ! never started
      if (self%count_rate == 0_int64) return          ! clock unavailable
      self%t_elapsed(i) = real(ticks - self%t_start(i), wp) &
           &              / real(self%count_rate, wp)
    end if

  end subroutine timer_measure


  subroutine timer_write(self, channel, headline)
    !***********************************************
    !* Print all finished timings to a Fortran     *
    !* I/O unit.                                   *
    !*                                             *
    !* Input:                                      *
    !*   channel  – Fortran unit number            *
    !*   headline – (optional) header line         *
    !***********************************************
    class(gfnff_timer), intent(in)           :: self
    integer,            intent(in)           :: channel
    character(len=*),   intent(in), optional :: headline

    integer  :: i
    real(wp) :: total

    if (present(headline)) then
      write(channel, '(/,1x,a)') trim(headline)
      write(channel, '(1x,a)')   repeat('-', 50)
    end if

    total = 0.0_wp
    do i = 1, self%n
      if (self%t_elapsed(i) < 0.0_wp) cycle       ! not finished
      total = total + self%t_elapsed(i)
      write(channel, '(1x,a,t36,f10.3,1x,a)') &
           & trim(self%labels(i)), self%t_elapsed(i), 's'
    end do

    if (self%n > 1) then
      write(channel, '(1x,a)') repeat('-', 50)
      write(channel, '(1x,a,t36,f10.3,1x,a)') 'total', total, 's'
    end if

  end subroutine timer_write

end module gfnff_timer_mod
