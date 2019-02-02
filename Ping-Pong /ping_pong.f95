program ping_pong
	
	! Include MPI module
	use mpi
	implicit none
	
	! Data declarations for MPI
	integer 	:: ierr ! error signal variable, Standard value - 0
	integer 	:: rank ! process ID (pid) / Number
	integer		:: nprocs ! number of processors
	
	! status variable - tells the status of send/ received calls
	! Needed for receive subroutine
	integer, dimension(MPI_STATUS_SIZE) :: status1
	
	! varaibles for neighbouring process
	integer :: left = 0, center = 0, right = 0
	
	! Test data: ping - data to send, pong - datya to receive
	integer		:: ping, pong
	
	! Initialize MPI
	call MPI_INIT(ierr)
	
	! Setup comminicator size
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
	
	! Setup rank/IDs for eac h process
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	
	! type main code here
	
	! setting up neighbouring processes
	left = rank -1
	center = rank
	right = rank +1
	
	! ping- pong game
	if (rank ==0) then
		! Set up a ping value and proceed
		ping = center +1
		
		!Snd/Pinga the data
		!Syntax: call MPI_SEND(start_address, count, datatype, destination pid/rank, tag, comminicator, ierr)
		call MPI_SEND(ping, 1, MPI_INTEGER, right, center, MPI_COMM_WORLD, ierr)
		print *, "Ping from rank", rank, "ping = ", ping
	
	
	else if (rank == nprocs-1) then
	
		! Receive/ Pong from process rank - 1
		! syntax call MPI_RECV(start_address, count, datatype, sourc, tag, communicator, status, ierr)
		call MPI_RECV(pong, 1, MPI_INTEGER, left, left, MPI_COMM_WORLD, status1, ierr)
		print *, "Pong from Rank", rank, "pong = ", pong
	
	else
	
		! Receive/ Pong from process rank - 1
		! syntax call MPI_RECV(start_address, count, datatype, sourc, tag, communicator, status, ierr)
		call MPI_RECV(pong, 1, MPI_INTEGER, left, left, MPI_COMM_WORLD, status1, ierr)
		print *, "Pong from Rank", rank, "pong = ", pong
		
		! modifying the data
		ping = pong +1
		
		!Snd/Pinga the data
		!Syntax: call MPI_SEND(start_address, count, datatype, destination pid/rank, tag, comminicator, ierr)
		call MPI_SEND(ping, 1, MPI_INTEGER, right, center, MPI_COMM_WORLD, ierr)
		print *, "Ping from rank", rank, "ping = ", ping

			
	end if
		
	! Finalize MPi
	call MPi_FINALIZE(ierr)
	
end program ping_pong
