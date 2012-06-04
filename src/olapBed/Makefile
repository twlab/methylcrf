

CC=g++

olapBed: olapBed.o trees/interval_tree.o
	$(CC) -o olapBed olapBed.o trees/interval_tree.o
f: olapBed
	./olapBed  s.bed t.bed
c: olapBed
	./olapBed -c s.bed t.bed
s: olapBed
	./olapBed -s s.bed t.bed
#	cat t.bed|./olapBed - 
clean: 
	-rm olapBed 2>/dev/null
	-rm *.o 2>/dev/null
