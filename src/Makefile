ODIR = ../obj/
BDIR = ../bin/
IDIR = ../img/
LIBS = -I/usr/include/eigen3 -I/usr/include/libcereal
CC = g++ -std=c++11 -pthread -Wall $(LIBS)
DEBUG = 
CFLAGS = -c -O2 -fPIC $(DEBUG)
LFLAGS = $(DEBUG)

OBJS = polymer.o stochastic.o spring.o angle.o \
       time_integration.o analysis.o rousetransform.o
OBJSFULL = $(patsubst %,$(ODIR)%,$(OBJS))

all: cpp python

# library installation on ubuntu
library_install:
	sudo apt-get install libeigen3-dev libcereal-dev python3-dev python3-sympy

# The C++ program 
cpp: $(BDIR)main $(BDIR)test 

$(BDIR)main : main.cpp main.h $(OBJSFULL) 
	$(CC) $(OBJSFULL) main.cpp -o $(BDIR)main $(LFLAGS) 
$(BDIR)test : unit_test.cpp $(OBJSFULL)
	$(CC) $(OBJSFULL) unit_test.cpp -o $(BDIR)test $(LFLAGS) 
$(ODIR)polymer.o : polymer.cpp
	$(CC) polymer.cpp -o $(ODIR)polymer.o $(CFLAGS)
$(ODIR)stochastic.o : stochastic.cpp 
	$(CC) stochastic.cpp -o $(ODIR)stochastic.o $(CFLAGS)
$(ODIR)analysis.o : analysis.cpp 
	$(CC) analysis.cpp -o $(ODIR)analysis.o $(CFLAGS)
$(ODIR)spring.o : spring.cpp 
	$(CC) spring.cpp -o $(ODIR)spring.o $(CFLAGS)
$(ODIR)angle.o : angle.cpp 
	$(CC) angle.cpp -o $(ODIR)angle.o $(CFLAGS)
$(ODIR)time_integration.o : time_integration.cpp 
	$(CC) time_integration.cpp -o $(ODIR)time_integration.o $(CFLAGS)
$(ODIR)rousetransform.o : rousetransform.cpp 
	$(CC) rousetransform.cpp -o $(ODIR)rousetransform.o $(CFLAGS)

# sudo locate Python.h to determine the right path
PLIB = -I/usr/include/python3.6m/ 
python: $(OBJSFULL) $(BDIR)_pypolymer.so

$(BDIR)_pypolymer.so: $(ODIR)pypolymer_wrap.o $(ODIR)main.o $(OBJSFULL)
	$(CC) $(LFLAGS) -shared $(ODIR)pypolymer_wrap.o $(ODIR)main.o \
	  $(OBJSFULL) -o $(BDIR)_pypolymer.so
	cp $(BDIR)pypolymer.py $(IDIR)pypolymer.py -rf
	cp $(BDIR)_pypolymer.so $(IDIR)_pypolymer.so -rf

$(ODIR)main.o: main.cpp main.h
	$(CC) $(CFLAGS) main.cpp -o $(ODIR)main.o

$(ODIR)pypolymer_wrap.o: pypolymer.i eigen.i analysis.i rousetransform.i
	rm $(BDIR)pypolymer.py || true
	swig -python -c++ $(LIBS) pypolymer.i
	mv pypolymer.py $(BDIR)pypolymer.py
	$(CC) $(CFLAGS) $(PLIB) pypolymer_wrap.cxx -o $(ODIR)pypolymer_wrap.o

clear: clean
clean:
	\rm $(BDIR)main $(BDIR)test $(BDIR)*.so || true
	\rm $(ODIR)*.o || true
	\rm *_wrap.cxx || true
