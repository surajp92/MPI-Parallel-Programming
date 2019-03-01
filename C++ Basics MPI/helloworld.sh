clear
ls
# Remove old executables
rm *.exe
ls

# COmpile and Build the executable
/usr/bin/mpic++ -o helloworld.exe helloworld.cpp
ls
# Run the executable
/usr/bin/mpirun -n 4 ./helloworld.exe
