CILK?=0
NATIVE?=1
OPT?=3
SANITIZE?=0

CFLAGS := -Wall -Wextra -O$(OPT) -g  -std=c++20 -gdwarf-4 -fno-exceptions -Wno-unknown-pragmas -Wno-comment

ifeq ($(NATIVE),1)
CFLAGS += -march=native
endif

ifeq ($(CILK),1)
CFLAGS += -fopencilk
endif

ifeq ($(SANITIZE),1)
ifeq ($(CILK),1)
CFLAGS += -fsanitize=cilk,undefined,address -fno-omit-frame-pointer
else
CFLAGS += -fsanitize=undefined,address -fno-omit-frame-pointer
endif
endif

DEFINES := -DCILK=$(CILK)

all: parspmv both_d spmm_dall spmm_a spmm_sall


seqsym: sym_spmv_test.cpp csbsym.cpp csbsym.h utility.h friends.h SSEspmv.o
	$(CXX) $(CFLAGS) $(DEFINES) -o seqsym sym_spmv_test.cpp SSEspmv.o

parsym: sym_spmv_test.cpp csbsym.cpp csbsym.h utility.h friends.h SSEspmv.o
	$(CXX) $(CFLAGS) $(DEFINES) -o parsym sym_spmv_test.cpp SSEspmv.o 

symanal: sym_spmv_test.cpp csbsym.cpp csbsym.h utility.h friends.h SSEspmv.o
	$(CXX) $(CFLAGS) $(DEFINES) -o symanal sym_spmv_test.cpp SSEspmv.o

seqspmv: csb_spmv_test.cpp bicsb.cpp bicsb.h bmcsb.cpp bmcsb.h friends.h utility.h SSEspmv.o
	$(CXX) $(CFLAGS) $(DEFINES) -o seqspmv csb_spmv_test.cpp SSEspmv.o

parspmv: csb_spmv_test.cpp bicsb.cpp bicsb.h bmcsb.cpp bmcsb.h friends.h utility.h SSEspmv.o 
	$(CXX) $(CFLAGS) $(DEFINES) -o parspmv csb_spmv_test.cpp SSEspmv.o

parspmv_nobm: csb_spmv_test.cpp bicsb.cpp bicsb.h friends.h utility.h
	$(CXX) $(CFLAGS) $(DEFINES) -DNOBM -o parspmv_nobm csb_spmv_test.cpp

parspmvt: csb_spmvt_test.cpp bicsb.cpp bicsb.h utility.h friends.h
	$(CXX) $(CFLAGS) $(DEFINES) -o parspmvt csb_spmvt_test.cpp

both_d:	both_test.cpp bicsb.cpp bicsb.h utility.h friends.h
	$(CXX) $(CFLAGS) $(DEFINES) -o both_d both_test.cpp

both_s:	both_test.cpp bicsb.cpp bicsb.h utility.h friends.h
	$(CXX) $(CFLAGS) $(DEFINES) -DSINGLEPRECISION -o both_s both_test.cpp

spmm_dall:	spmm_test.cpp bicsb.cpp bicsb.h utility.h friends.h
	for number in 4 8 12 16 24 32 40 48 56 64; do \
		echo "$(CXX) $(CFLAGS) $(DEFINES) -DRHSDIM=$$number -o spmm_d$$number spmm_test.cpp"; \
		$(CXX) $(CFLAGS) $(DEFINES) -DRHSDIM=$$number -o spmm_d$$number spmm_test.cpp; \
	done;

spmm_a:	spmm_test.cpp bicsb.cpp bicsb.h utility.h friends.h
	$(CXX) $(CFLAGS) $(DEFINES) -DSINGLEPRECISION -S -fcode-asm -vec_report6 spmm_test.cpp

spmm_sall:	spmm_test.cpp bicsb.cpp bicsb.h utility.h friends.h
	for number in 4 8 12 16 24 32 40 48 56 64; do \
		echo "$(CXX) $(CFLAGS) $(DEFINES) -DSINGLEPRECISION -DRHSDIM=$$number -o spmm_s$$number spmm_test.cpp"; \
		$(CXX) $(CFLAGS) $(DEFINES) -DSINGLEPRECISION -DRHSDIM=$$number -o spmm_s$$number spmm_test.cpp; \
	done;

SSEspmv.o: SSEspmv.cpp
	$(CXX) $(CFLAGS) $(DEFINES) -c SSEspmv.cpp	

clean:	
	rm -f seqspmv
	rm -f seqsym
	rm -f parspmv
	rm -f parsym 
	rm -f parspmvt
	rm -f parspmv_nobm
	for number in 8 16 24 32 40 48 56 64; do \
		rm -f spmm_s$$number;\
		rm -f spmm_d$$number;\
	done;
	rm -f *.o
