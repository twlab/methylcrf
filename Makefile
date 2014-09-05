
# g++ will support -std=c++11
CC=g++
CXXFLAGS= -W -Wall -Wextra -pedantic -std=c++0x

olapBed: olapBed.o trees/interval_tree.o
	$(CC) $(CXXFLAGSS) -o olapBed olapBed.o trees/interval_tree.o

f: olapBed
	./olapBed  s.bed t.bed
c: olapBed
	./olapBed -c s.bed t.bed
s: olapBed
	./olapBed -s s.bed t.bed
#	cat t.bed|./olapBed - 
