# @Author: Suraj Pawar <user1>
# @Date:   2019-02-09T14:32:19-06:00
# @Email:  supawar@okstate.edu
# @Filename: 2dLaplace.sh
# @Last modified by:   user1
# @Last modified time: 2019-02-09T14:32:24-06:00



clear
ls
# Remove old executables
rm *.exe
ls

# COmpile and Build the executable
/usr/bin/mpif90 -fbounds-check -o 2dLaplace.exe 2dLaplace.f95
ls
# Run the executable
/usr/bin/mpirun -n 6 ./2dLaplace.exe
