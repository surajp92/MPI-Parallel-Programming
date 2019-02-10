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
real*8, dimension(:), allocatable	:: u_final_local
real*8, parameter					:: x_left = -2.0, x_right = 2.0
integer, parameter					:: n_global = 321
integer								:: n_local
real*8								:: dx, x
integer								:: i, k
integer, parameter					:: istart = 1, kstart = 1
integer								:: i_global_low, i_global_high


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

! Call update subroutine to calculate velocity field
call update(rank, nprocs, n_local, n_global, u_initial_local, dx, x_left, dt, T)

! Call collect subroutine to collect velocity from all processors and store it in 0th processor
call collect(rank, nprocs, n_local, n_global, u_initial_local, dx, x_left)

! Finalize MPI
call MPI_FINALIZE(ierr)

deallocate(u_initial_local)

end program inviscid_burgers

subroutine update(rank, nprocs, n_local, n_global, u_initial_local, dx, x_left, dt, T)

use mpi

implicit none

integer								:: i_local_low, i_local_high
integer								:: i_global_low, i_global_high
integer								:: i_local, i_global
integer								:: n_local, n_global
integer								:: k
real*8								:: dx, x, x_left
real*8								:: u_initial_local(n_local)
real*8								:: u_final_local(n_local)
real*8								:: dt, T, lambda
integer								:: time_steps

! Data declarations for MPI
integer 	:: ierr ! error signal variable, Standard value - 0
integer 	:: rank ! process ID (pid) / Number
integer		:: nprocs ! number of processors

! status variable - tells the status of send/ received calls
! Needed for receive subroutine
integer, dimension(MPI_STATUS_SIZE) :: status

! tags for send/ receive
integer, parameter					:: right = 10
integer, parameter					:: left = 20

time_steps = T/dt
lambda = dt/dx

i_global_low  = (rank		*(n_global-1))/nprocs
i_global_high = ((rank+1)	*(n_global-1))/nprocs

if (rank > 0) then
	i_global_low = i_global_low - 1
end if

i_local_low = 0
i_local_high = i_global_high - i_global_low

! Assign initial condition
do i_global = i_global_low, i_global_high
	x = x_left + dx*(i_global)
	i_local = i_global - i_global_low
	if (x .le. -1 .or. x.ge.1) then
		u_initial_local(i_local + 1) = 0.0
	else
		u_initial_local(i_local + 1) = 1-abs(x)
	end if
end do

do k = 2,time_steps

print *, rank, k

!###############  calculate results at new time step ###############!
	do i_local = i_local_low + 1 , i_local_high - 1
		u_final_local(i_local+1) = u_initial_local(i_local+1) &
								  -0.5*lambda*(u_initial_local(i_local+1+1)-u_initial_local(i_local+1-1))
	end do


!############### communicate to share common points between two processors ###############!

! Comminication to left processor and from right processor
	if (rank > 0) then
		!Send (i_local_low+2) to (i_local_high+1) of left
		!Syntax: call MPI_SEND(start_address, count, datatype, destination pid/rank, tag, comminicator, ierr)
		call MPI_SEND(u_final_local(i_local_low+2), 1, MPI_DOUBLE_PRECISION, rank-1, left, MPI_COMM_WORLD, ierr)

		! receive (i_local_high)
		! syntax call MPI_RECV(start_address, count, datatype, sourc, tag, communicator, status, ierr)
	 	call MPI_RECV(u_final_local(i_local_low+1), 1, MPI_DOUBLE_PRECISION, rank-1, right, MPI_COMM_WORLD, status, ierr)
	else
		! Boundary condiiton on left side for processor 0
		u_final_local(i_local_low+1) = 0.0
	end if

! Communication from left processor and to the right processor
	if (rank < nprocs-1) then
		! receive (i_local_low+2) of left to (i_local_high+1)
		! syntax call MPI_RECV(start_address, count, datatype, sourc, tag, communicator, status, ierr)
	 	call MPI_RECV(u_final_local(i_local_high+1), 1, MPI_DOUBLE_PRECISION, rank+1, left, MPI_COMM_WORLD, status, ierr)

	 	!Send (i_local_high) of left to (i_local_low+1) of right
		!Syntax: call MPI_SEND(start_address, count, datatype, destination pid/rank, tag, comminicator, ierr)
		call MPI_SEND(u_final_local(i_local_high), 1, MPI_DOUBLE_PRECISION, rank+1, right, MPI_COMM_WORLD, ierr)
	else
		! Boundary condiiton on right side for nth processor
		u_final_local(i_local_high+1) = u_final_local(i_local_high)
	end if

!###############  update initial field with final field	###############!
	do i_local = i_local_low, i_local_high
		u_initial_local(i_local+1) = u_final_local(i_local+1)

	end do
end do


!print *, 'Rank', rank
!print *, i_global_low, i_global_high
!print *, i_local_low, i_local_high
!print *, u_initial_local
!1 format(f2.2)
!print (*, Format), 'Initial velocity field', u_initial_local

return
end subroutine update

subroutine collect(rank, nprocs, n_local, n_global, u_initial_local, dx, x_left)

use mpi

implicit none

integer								:: i_local_low, i_local_high
integer								:: i_global_low, i_global_high
integer								:: i_local, i_global
integer								:: n_local, n_global
real*8								:: u_initial_local(n_local)
real*8, dimension(:), allocatable	:: u_global, x_global
real*8, dimension(:), allocatable	:: u_inital_global
integer								:: procs, i
integer								:: n_local_procs
real*8								:: dx, x_left

! Data declarations for MPI
integer 	:: ierr ! error signal variable, Standard value - 0
integer 	:: rank ! process ID (pid) / Number
integer		:: nprocs ! number of processors

! MPI send/ receive arguments
integer								:: buffer(2)
integer, parameter					:: collect1 = 10
integer, parameter					:: collect2 = 20


! status variable - tells the status of send/ received calls
! Needed for receive subroutine
integer, dimension(MPI_STATUS_SIZE) :: status

i_global_low  = (rank		*(n_global-1))/nprocs
i_global_high = ((rank+1)	*(n_global-1))/nprocs

if (rank > 0) then
	i_global_low = i_global_low - 1
end if

i_local_low = 0
i_local_high = i_global_high - i_global_low

	allocate(x_global(1:n_global))
	allocate(u_inital_global(1:n_global))

	do i_global = 1, n_global
		x_global(i_global) = x_left + dx*(i_global-1)
		if (x_global(i_global) .le. -1 .or. x_global(i_global) .ge. 1) then
			u_inital_global(i_global) = 0.0
		else
			u_inital_global(i_global) = 1-abs(x_global(i_global))
		end if
	end do

if (rank == 0) then
	allocate(u_global(1:n_global))


	do i_local = i_local_low, i_local_high
		i_global = i_global_low + i_local - i_local_low
		u_global(i_global+1) = u_initial_local(i_local+1)
	end do

	do procs = 1,nprocs-1

		! receive global_low and number_local points in the processor
		! syntax call MPI_RECV(start_address, count, datatype, sourc, tag, communicator, status, ierr)
		call MPI_RECV(buffer, 2, MPI_INTEGER, procs, collect1, MPI_COMM_WORLD, status, ierr)

		i_global_low = buffer(1)
		n_local_procs = buffer(2)

		if ( i_global_low < 0 ) then
			write ( *, '(a,i6)' ) '  Illegal I_GLOBAL_LO = ', i_global_low
			call MPI_Finalize ( ierr )
			stop 1
		else if ( n_global <= i_global_low + n_local_procs - 1 ) then
        write ( *, '(a,i6)' ) '  Illegal I_GLOBAL_LO + N_LOCAL = ', &
			i_global_low + n_local_procs
			call MPI_Finalize ( ierr )
			stop 1
		end if

		! receive u_local from the processor and add it to u_global
		! syntax call MPI_RECV(start_address, count, datatype, sourc, tag, communicator, status, ierr)
	 	call MPI_RECV(u_global(i_global_low+1), n_local_procs, MPI_DOUBLE_PRECISION, procs, collect2, MPI_COMM_WORLD, status, ierr)

	end do

	open(unit = 4, file = 'solution.dat')
	write(4,*), 'Varaibles = "X""Initial""Final"'
	do i = 1, n_global
		write(4,*), x_global(i), u_inital_global(i), u_global(i)
	end do

	! print *, u_global
	! print *, i_global_low, n_local_procs
	! Writing the final solution
	deallocate(u_global)

else

	buffer(1) = i_global_low
	buffer(2) = n_local

	!Send global_low and n_local for the processor
	!Syntax: call MPI_SEND(start_address, count, datatype, destination pid/rank, tag, comminicator, ierr)
	call MPI_SEND(buffer, 2, MPI_INTEGER, 0, collect1, MPI_COMM_WORLD, ierr)

	!Send u_local from the processor and add it to u_global
	!Syntax: call MPI_SEND(start_address, count, datatype, destination pid/rank, tag, comminicator, ierr)
	call MPI_SEND(u_initial_local, n_local, MPI_DOUBLE_PRECISION, 0, collect2, MPI_COMM_WORLD, ierr)

end if

return
end subroutine collect
