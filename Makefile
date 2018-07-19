


all: methylCRF.pl

methylCRF.pl: MRE_norm.pl medip_norm.sh crfasgd bed2avgwinbin.sh dir2frag.sh sam2bed.pl
	ln -s src/methylCRF.pl ./

MRE_norm.pl:
	ln -s src/MRE_norm.pl ./

bed2avgwinbin.sh:
	ln -s src/bed2avgwinbin.sh ./

dir2frag.sh:
	ln -s src/dir2frag.sh ./

medip_norm.sh:
	ln -s src/medip_norm.sh ./

sam2bed.pl:
	ln -s src/sam2bed.pl ./

crfasgd: 
	cd src/sgd/crf; make; cd ../../../;
	ln -s src/sgd/crf/crfasgd ./

# mapBed: bedtools should be installed
#	ln -s /opt/apps/bedtools/2.25.0/bin/mapBed ./


clean:
	-rm methylCRF.pl MRE_norm.pl sam2bed.pl medip_norm.sh crfasgd bed2avgwinbin.sh dir2frag.sh 2>/dev/null
	cd src/sgd/crf; make clean; cd ../../../;


.PHONY: all clean


