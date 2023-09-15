rm *.o
opt='-Wall -g'
CPP=g++
exe=d2ts
$CPP -c d2ts.cpp tabu_search.cpp reads.cpp main.cpp $opt 
$CPP -o $exe d2ts.o tabu_search.o reads.o main.o $opt
