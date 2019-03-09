ROOT_LIBS=$(shell root-config --libs)
ROOT_INCLUDE=$(shell root-config --cflags)

armaLib=../../../arma/install/lib64/
armaInc=../../../arma/install/include/

readcov.o: readcov.cpp
	g++ -c $(ROOT_INCLUDE) -I$(armaInc)  $^  $(ROOT_LIBS) -L$(armaLib)  -larmadillo  -o $@    

readcov: readcov.cpp
	g++  $(ROOT_INCLUDE) -I$(armaInc)  $^  $(ROOT_LIBS) -L$(armaLib) -larmadillo  -o $@    

