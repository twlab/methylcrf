


all: methylCRF.pl

methylCRF.pl: medip_norm.sh crfsgd bed2avgwinbin.sh dir2frag.sh
	ln -s src/methylCRF.pl ./

medip_norm.sh:
	ln -s src/medip_norm.sh ./

crfsgd: 
	cd src/sgd/crf; make; cd ../../../;
	ln -s src/sgd/crf/crfsgd ./

bed2avgwinbin.sh: olapBed
	ln -s src/bed2avgwinbin.sh ./

dir2frag.sh: olapBed
	ln -s src/dir2frag.sh ./

olapBed: 
	cd src/olapBed; make; cd ../../;
	ln -s src/olapBed/olapBed ./

clean:
	-rm methylCRF.pl medip_norm.sh crfsgd bed2avgwinbin.sh dir2frag.sh olapBed 2>/dev/null
	cd src/sgd; make clean; cd ../../;
	cd src/olapBed; make clean; cd ../


.PHONY: all clean



