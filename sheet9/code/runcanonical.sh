g++ -o main main.cpp -fopenmp

mode=1 # mode 1 is canonical, mode 0 in microcanonical
N1d=8 # number of particles per dimension
Nsteps=60000 # number of simulation steps
T1=80 # temperature simulation 1
T2=70 # temperature simulation 2
T3=395 # temperature simulation 3
T4=405 # temperature simulation 4

./main $mode $T1 $N1d $Nsteps > "temp$T1.txt"
./main $mode $T2 $N1d $Nsteps > "temp$T2.txt"
./main $mode $T3 $N1d $Nsteps > "temp$T3.txt"
./main $mode $T4 $N1d $Nsteps > "temp$T4.txt"

python3 plot.py 1