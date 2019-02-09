clear
ls
# Remove old executables
rm *.exe
ls

# COmpile and Build the executable
/usr/bin/mpif90 -fbounds-check -o mpi_inviscid_burgers.exe inviscid_burgers.f95
ls
# Run the executable
/usr/bin/mpirun -n 4 ./mpi_inviscid_burgers.exe
