


all: methylCRF.pl

methylCRF.pl: MRE_norm.pl medip_norm.sh crfasgd bed2avgwinbin.sh dir2frag.sh sam2bed.pl
	ln -s src/methylCRF.pl ./

MRE_norm.pl:
	ln -s src/MRE_norm.pl ./

bed2avgwinbin.sh: mapBed
	ln -s src/bed2avgwinbin.sh ./

dir2frag.sh: mapBed
	ln -s src/dir2frag.sh ./

medip_norm.sh: mapBed
	ln -s src/medip_norm.sh ./

sam2bed.pl:
	ln -s src/sam2bed.pl ./

crfasgd: 
	cd src/sgd/crf; make; cd ../../../;
	ln -s src/sgd/crf/crfasgd ./

mapBed:
	ln -s /usr/local/bin/mapBed ./

# olapBed: 
# 	cd src/olapBed; make; cd ../../;
# 	ln -s src/olapBed/olapBed ./

clean:
	-rm methylCRF.pl MRE_norm.pl sam2bed.pl medip_norm.sh crfasgd bed2avgwinbin.sh dir2frag.sh mapBed 2>/dev/null
	cd src/sgd; make clean; cd ../../;
	# cd src/olapBed; make clean; cd ../


.PHONY: all clean


