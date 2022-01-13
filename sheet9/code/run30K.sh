g++ -o main main.cpp -fopenmp

mode=1 # mode 1 is canonical, mode 0 in microcanonical
N1d=8 # number of particles per dimension
Nsteps=60000 # number of simulation steps
T1=30 # temperature simulation 1

./main $mode $T1 $N1d $Nsteps > "temp$T1.txt"

python3 plot.py 1