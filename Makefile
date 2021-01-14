GCCOPT = -O2 -fno-rtti -fno-exceptions # -ftree-vectorize
INTELOPT = -O2 -no-ipo -fno-rtti -fno-exceptions -parallel -restrict -std=c++11 -xAVX -no-prec-div #-fno-inline-functions
DEB = -g -DNOBM -O0 -parallel -restrict -std=c++11 


seqspmv: csb_spmv_test.cpp bicsb.cpp bicsb.h friends.h utility.h 
	g++ $(INCADD) $(GCCOPT) -o seqspmv csb_spmv_test.cpp 


clean:	
	rm -f seqspmv
	rm -f *.o
