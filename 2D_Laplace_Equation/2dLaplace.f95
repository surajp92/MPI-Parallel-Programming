!------------------------------------------------------------------
! @Author: Suraj Pawar <user1>
! @Date:   2019-02-09T14:27:54-06:00
! @Email:  supawar@okstate.edu
! @Last modified by:   user1
! @Last modified time: 2019-02-09T14:27:54-06:00
!------------------------------------------------------------------

!------------------------------------------------------------------
! This program solves Two-dimensional Laplace equation in square
! domain using MPI.
! To compile the code you can use the bash script available along
! with the code. You can change the number of cores in .sh file
! Use bash 2dLaplace.sh to run the code
!------------------------------------------------------------------

program laplace

implicit none

include 'mpif.h'

! Global domain parameters
integer, parameter          ::           nx_global = 256
integer, parameter          ::           ny_global = 256
real*8, parameter           :: max_error_tolerance = 1.0e-6
real*8, parameter           ::      max_iterations = 50000
real*8                      ::    max_error_global = 100.0

! Local parameters for processors
integer                     :: nx, ny
real*8                      :: start_time, stop_time
integer                     :: i,j
integer                     :: iteration_count
real*8                      :: max_error_local
real*8                      :: h, xp, yp, yp_loc_min

! Temperature Solution and x, y location for each grid
real*8, allocatable, dimension(:,:)     :: temperature
real*8, allocatable, dimension(:,:)     :: temperature_old
real*8, allocatable, dimension(:,:)     :: temperature_exact
real*8, allocatable, dimension(:,:)     :: x, y

! MPI reated variables
integer 	                           :: ierr ! error signal variable, Standard value - 0
integer 	                           :: rank ! process ID (pid) / Number
integer		                           :: nprocs ! number of processors
integer, dimension(MPI_STATUS_SIZE)    :: status ! status for receive
integer, parameter					   :: id_top_to_bottom = 10 ! tag for send/ receive
integer, parameter					   :: id_bottom_to_top = 20 ! tag for send/ receive
integer, dimension(4)                  :: req = MPI_REQUEST_NULL
integer, dimension(MPI_STATUS_SIZE, 4) :: status_array

!------------------------------------------------------------------
! MPI Starting

call MPI_INIT(ierr) ! Initialize MPI

call MPI_COMM_SIZE(MPI_COMM_WORLD, &
                           nprocs, &    ! number of processors
                           ierr     )   ! ierr = 0 if successful

call MPI_COMM_RANK(MPI_COMM_WORLD, &
                             rank, &    ! rank of the processor
                             ierr   )   ! ierr = 0 if successful

if (rank == 0) then
    write(*,*) 'No. of processor = ', nprocs, 'Rank = ', rank
    write(*,*)
    write(*,*) 'Nx_Global = ', nx_global
    write(*,*) 'Ny_Global = ', ny_global
    write(*,*) 'Total size = ', ny_global*nx_global
end if

!------------------------------------------------------------------
! Allocate x and y position to each grid point.
! x = 0 corresponds to i =0 and x = 1 corresponds to i = nx_global + 1
! The domain is partitioned in y direction and each domain is assigned
! to individual processor
! Each domain has ghost nodes on top and bottom side. So, they have
! (ny+2) points in y- direction
! (J = 0)_(rank) = (J = ny+1)_(rank-1)
! (J = ny+1)_(rank) = (J = 0)_(rank+1)

! local size
nx = nx_global
ny = ny_global/nprocs
!print *, nx, ny+1

! grid spacing size
h = 1.0/(ny_global)

! allocate the local array for temperature grid positions
allocate(      temperature(0:nx+1, 0:ny+1))
allocate(  temperature_old(0:nx+1, 0:ny+1))
allocate(temperature_exact(0:nx+1, 0:ny+1))
allocate(                x(0:nx+1, 0:ny+1))
allocate(                y(0:nx+1, 0:ny+1))

! grid position data
do j = 0, ny+1

    ! y-coordinate
    yp_loc_min = rank*(1.0/nprocs) ! dimesnion in y-direction is 1.0
            yp = yp_loc_min + h*(j-1) ! (j-1) = bottom ghost, j = 0 is ghost node at -h

    do i = 0, nx+1
        x(i,j) = real(i)/real(nx+1) ! x(0) = 0.0, x(nx+1) = 1.0
        y(i,j) = yp
    end do

end do

!------------------------------------------------------------------
! Calculate initial, boundary condition and exact condition

call initialize (rank, nprocs, x, y, nx, ny, temperature_old, temperature_exact)
temperature = temperature_old

!------------------------------------------------------------------
! start recording actual execution time

call cpu_time(start_time)

! initialize the iteration count
iteration_count = 1

! continue iteration until max_error_global > = max_error_tolerance or
! iteration_count <= max_iterations

do while (max_error_global >= max_error_tolerance .and. iteration_count <= max_iterations)

! compute solution at next time step using explicit scheme
    do j = 1, ny
        ! skip bottom boundary at j = 1
        if (rank == 0 .and. j ==1) cycle
        ! do calculation for all other points in local domain
        ! skip left and right boundary
        do i = 1, nx
            temperature(i,j) = 0.25*(temperature_old(i-1,j) &
                                   + temperature_old(i+1,j) &
                                   + temperature_old(i,j+1) &
                                   + temperature_old(i,j-1))
        end do

    end do


!------------------------------------------------------------------
! communicate the solution between ghost nodes

    ! send solution of top boundary of local domain (rank 0,1..) (j = ny) to
    ! bottom ghost nodes of upper domain ((rank+1) 1,2..) (j = 0)
    if (rank < nprocs-1) then
        call MPI_Isend(temperature(1,ny), & ! starting address of data to be send
                                      nx, & ! number of bytes of data to be send
                    MPI_DOUBLE_PRECISION, & ! datatype of data to be send
                                 rank+1 , & ! target processor rank
                        id_top_to_bottom, & ! tag of the send signal
                          MPI_COMM_WORLD, & !
                                  req(1), & !
                                    ierr)   ! ierr = 0 if successful
    end if

    ! send solution of bottom boundary of local domain (rank 1,2..) (j = 1) to
    ! top ghost nodes of lower domain ((rank-1) 0,1..) (j = ny+1)
    if (rank /= 0) then
        call MPI_Isend(temperature(1,1), & ! starting address of data to be send
                                     nx, & ! number of bytes of data to be send
                   MPI_DOUBLE_PRECISION, & ! datatype of data to be send
                                rank-1 , & ! target processor rank
                       id_bottom_to_top, & ! tag of the send signal
                         MPI_COMM_WORLD, & !
                                 req(2), & !
                                   ierr)   ! ierr = 0 if successful
    end if

    ! receive data from bottom domain (rank-1 0,1) (j = ny)  and
    ! store it in bottom ghost nodes (j = 0)
    if (rank /= 0) then
        call MPI_Irecv(temperature_old(1,0), & ! starting address of data to be receive
                                     nx, & ! number of bytes of data to be receive
                   MPI_DOUBLE_PRECISION, & ! datatype of data to be receive
                                rank-1 , & ! source processor rank
                       id_top_to_bottom, & ! tag of the send signal
                         MPI_COMM_WORLD, & !
                                 req(3), & !
                                   ierr)   ! ierr = 0 if successful
    end if

    ! receive data from tom domain (rank+1 1,2) (j = 1)  and
    ! store it in top ghost nodes (j = ny+1)
    if (rank /= nprocs-1) then
        call MPI_Irecv(temperature_old(1,ny+1), & ! starting address of data to be send
                                        nx, & ! number of bytes of data to be send
                      MPI_DOUBLE_PRECISION, & ! datatype of data to be send
                                   rank+1 , & ! source processor rank
                          id_bottom_to_top, & ! tag of the send signal
                            MPI_COMM_WORLD, & !
                                    req(3), & !
                                      ierr)   ! ierr = 0 if successful
    end if


    ! wait till Isend and Irecv above is completed
    call MPI_Waitall(4, req, status_array, ierr)

    ! compute maximum local error between temperature and temperature_old
    max_error_local = 0.0

    do j = 1, ny
        if (rank == 0 .and. j ==1) cycle
        do i = 1,nx
            max_error_local = max(abs(temperature(i,j)-temperature_old(i,j)), max_error_local)
            temperature_old(i,j) = temperature(i,j)
        end do
    end do

    !print *, 'Final initial condition'
    !do j = 0, ny+1
    !    print *, rank, 'new', j, temperature(:,j)
        ! print *, rank, 'new', j, temperature(:,j)
    !end do

    ! reduce maximum_error_local in each domain to max_error_global
    ! and store it in rank = 0 processor

    call MPI_Reduce(max_error_local, & ! variable to be collected from all processors
                   max_error_global, & ! variable in which result is to be stored
                                  1, & ! number of values
               MPI_DOUBLE_PRECISION, & ! datatype of variable
                            MPI_MAX, & ! type of operation
                                  0, & ! rank of processor in which variable is to be stored
                     MPI_COMM_WORLD, & !
                                ierr)  ! ierr = 0 if successful

    ! Broadcast max_error_global to all processors

    call MPI_Bcast(max_error_global, & ! variable to be broadcasted
                                  1, & ! number of values
               MPI_DOUBLE_PRECISION, & ! datatype of variabel
                                  0, & ! braodcast from procs with rank = 0
                     MPI_COMM_WORLD, & !
                                ierr)  ! ierr = 0 if successful

    ! display iteration count and max_error_global at every iteration_count
    if (rank == 0 .and. mod(iteration_count,10) ==0) then
        write(*,*), 'Iter = ', iteration_count, 'Maximum Error = ', max_error_global
    end if

    iteration_count = iteration_count + 1
end do

! wait till all processors come to this points
! this is for accurate timing and clean output

call MPI_Barrier(MPI_COMM_WORLD, & !
                            ierr) ! ierr = 0 if successful

! store stop time

call cpu_time(stop_time)

! print total time, total iterations, and maximum global Error
if (rank == 0) then
    print*, 'Total Iteration Count = ', iteration_count, 'Maxium Error = ', max_error_global
    print*, 'Total Time = ', stop_time-start_time, 'seconds', 'No. of nodes', nx*ny
end if

write(rank+1000,*) '  Total time = ', stop_time-start_time, ' seconds.', ' nnodes=',nx*ny

temperature = temperature_old
call write_tecplot_file(x,y,temperature,temperature_exact,nx,ny,rank,nprocs)

! MPI final calls
call MPI_FINALIZE(ierr) ! ierr = 0 if successful

end program laplace

!------------------------------------------------------------------
! initialize subroutine calculates initial condition, boundary
! consition and exact Solution
! input: x(:,:), y(:,:) = x, y coordinate of each grid points
!                nx, ny = size of local domain in x and y direction
!                  Rank = rak of the processor
!                nprocs = number of processors
! output:   temperature_old = initial temperature field\
!         temperature_exact = exact solution
!------------------------------------------------------------------

subroutine initialize(rank, nprocs, x, y, nx, ny, temperature_old, temperature_exact)

implicit none

integer, intent(in)                             :: nx, ny, rank, nprocs
real*8, intent(in), dimension(0:nx+1, 0:ny+1)   :: x, y
real*8, intent(out), dimension(0:nx+1, 0:ny+1)  :: temperature_old
real*8, intent(out), dimension(0:nx+1, 0:ny+1)  :: temperature_exact

integer                                         :: i, j
real*8                                          :: xp, yp
real*8, parameter                               :: pi = 3.141592653589793238

! calculate exact solutuion and assign it to initial condition
do j = 0, ny+1

    do i = 0, nx+1
        xp = x(i,j)
        yp = y(i,j)
        temperature_exact(i,j) = (sinh(pi*xp)*sin(pi*yp) + sinh(pi*yp)*sin(pi*xp))/sinh(pi)
        temperature_old(i,j) = temperature_exact(i,j)
    end do

end do

! change initial condiiton for interior points. Do not change at boundaries
! exact solution at boundaries is the boundary condition


do j = 0, ny+1
    ! boundary condition on the bototm side of the domain
    ! j = 0 is ghost and j = 1 is lower boundarey
    if (rank == 0 .and. j<2) cycle

    ! boundary condition on the top side of the domain
    ! j = ny+1 is upper boundarey
    if (rank == nprocs-1 .and. j == ny+1) cycle

    ! skip i = 0 (left boundary) and i = nx+1 (right boundary)
    do i = 1, nx
        ! print *, rank, i, j
        temperature_old(i,j) = 0.0
    end do

end do

end subroutine initialize

!------------------------------------------------------------------
! subroutine for writing results into tecplot which can be plotted

subroutine write_tecplot_file(x,y,temperature,temperature_exact,nx,ny,rank,nprocs)

implicit none

integer , intent(in)                           :: nx, ny, rank,nprocs
real*8, intent(in), dimension(0:nx+1,0:ny+1) ::  x,  y
real*8, intent(in), dimension(0:nx+1,0:ny+1) :: temperature,temperature_exact

integer                                        :: i, j
character(80)                                  :: char_temp, filename

!---------------------------------------------------------------------------
! Store the value of 'rank' as a character in the character variable 'char_temp'.
 write(char_temp,'(i5)') rank

!---------------------------------------------------------------------------
! Define the file name
 filename = "tecplot_" // trim(adjustl(char_temp)) // '.dat'

!---------------------------------------------------------------------------
! Open the file and start writing the data!
 open(unit=1, file=filename, status="unknown")

!---------------------------------------------------------------------------
! Header
 write(1,*) 'title =', '"Partition_'// trim(adjustl(char_temp)), '"'
 write(1,*) 'variables = "x", "y", "Temperature", "Temperature exact" "T-Texact"'

 if (rank==0) then
  write(1,*) 'zone T=','Partition_'// trim(adjustl(char_temp)),' i=', nx+2, 'j=', ny+1, 'f=point'
 else
  write(1,*) 'zone T=','Partition_'// trim(adjustl(char_temp)),' i=', nx+2, 'j=', ny+2, 'f=point'
 endif

!---------------------------------------------------------------------------
! Note: Exclude the bottom ghost nodes in rank=0, which lie outside the domain, and so
!       are not used in the iteration and completely irrelevant to the soltuion.
! Write grid data (including the ghost points at j = 0 and ny+1)

 do j = 0, ny+1

  if (rank==0 .and. j==0) cycle

  do i = 0, nx+1

   write(1,*) x(i,j), y(i,j), temperature(i,j), temperature_exact(i,j), &
              temperature(i,j) - temperature_exact(i,j) ! <- error

  end do

 end do

!---------------------------------------------------------------------------
! Close the file
 close(1)

end subroutine write_tecplot_file
