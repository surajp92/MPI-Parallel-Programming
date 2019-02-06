clear
ls
# Remove old executables
rm *.exe
ls

# COmpile and Build the executable
/usr/bin/mpif90 -g -o mpi_inviscid_burgers inviscid_burgers.f95
ls
# Run the executable
/usr/bin/mpirun -np 2 xterm -e gdb mpi_inviscid_burgers
