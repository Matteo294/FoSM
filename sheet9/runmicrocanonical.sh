g++ -o main main.cpp -fopenmp

mode=0 # mode 1 is canonical, mode 0 in microcanonical
N1d=4 # number of particles per dimension
Nsteps=10000 # number of simulation steps
T=80 # temperature simulation (then translated into energy)

./main $mode $T $N1d $Nsteps > "micro_temp$T1.txt"

python3 plot.py 0