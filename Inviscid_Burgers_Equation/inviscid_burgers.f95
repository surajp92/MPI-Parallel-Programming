program inviscid_burgers

use mpi

implicit none

! Data declarations for MPI
integer 	:: ierr ! error signal variable, Standard value - 0
integer 	:: rank ! process ID (pid) / Number
integer		:: nprocs ! number of processors

! status variable - tells the status of send/ received calls
! Needed for receive subroutine
! integer, dimension(MPI_STATUS_SIZE) :: status1

real*8, parameter					:: dt=0.00125, T = 0.25
real*8, dimension(:), allocatable	:: u_initial_local !, u_next
real*8, parameter					:: x_left = -2.0, x_right = 2.0
integer*8, parameter				:: n_global = 9
integer*8							:: n_local
real*8								:: dx, x
integer*8							:: i, k
integer*8, parameter				:: istart = 1, kstart = 1
integer*8							:: i_global_low, i_global_high


dx = (x_right-x_left)/(n_global-1)

! Initialize MPI
call MPI_INIT(ierr)
	
! Setup comminicator size
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
	
! Setup rank/IDs for eac h process
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

i_global_low  = (rank		*(n_global-1))/nprocs
i_global_high = ((rank+1)	*(n_global-1))/nprocs

if (rank > 0) then
	i_global_low = i_global_low - 1
end if

n_local = i_global_high - i_global_low + 1
allocate(u_initial_local(1:n_local))  

call update(rank, nprocs, n_local, n_global, u_initial_local, dx, x_left)

!allocate(u_next(1:n_global))
!allocate(u_prev(1:n_global))

!do i = istart, n_global
!	x = x_left + (i-1)*dx
!	if (abs(x) .le. 1) then
!		u_prev(i) = 1-abs(x)
!	else
!		u_prev(i) = 0
!	end if	
!end do

!print *, u_prev

! Finalize MPi
call MPI_FINALIZE(ierr)
end program inviscid_burgers  

subroutine update(rank, nprocs, n_local, n_global, u_initial_local, dx, x_left)

use mpi

implicit none

integer*8							:: i_local_low, i_local_high
integer*8							:: i_global_low, i_global_high
integer*8							:: i_local, i_global
integer*8							:: n_local, n_global
integer								:: rank, nprocs
real*8								:: dx, x, x_left
real*8								:: u_initial_local(n_local)


i_global_low  = (rank		*(n_global-1))/nprocs
i_global_high = ((rank+1)	*(n_global-1))/nprocs

if (rank > 0) then
	i_global_low = i_global_low - 1
end if

i_local_low = 0
i_local_high = i_global_high - i_global_low

do i_global = i_global_low, i_global_high
	x = x_left + dx*(i_global)
	i_local = i_global - i_global_low
	u_initial_local(i_local + 1) = x
end do
print *, 'Rank', rank
print*, 'Initial velocity field', u_initial_local

end subroutine update 

