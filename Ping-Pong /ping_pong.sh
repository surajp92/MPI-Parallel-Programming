clear
ls
# Remove old executables
rm *.exe
ls

# COmpile and Build the executable
/usr/bin/mpif90 -o ping_pong.exe ping_pong.f95
ls
# Run the executable
/usr/bin/mpirun -n 16 ./ping_pong.exe
