CXX = g++-10
CC = gcc-10
VARNAME = value
CXXFLAGS = -Wall -g -O2 -fopenmp -std=c++17 -fPIC

main: main.o NumMethods.o quarkEOS.o nucEOS.o Conversions.o MCMC.o
	$(CXX) $(CXXFLAGS) -o main main.o NumMethods.o quarkEOS.o nucEOS.o Conversions.o MCMC.o
	
NumMethods.o: NumMethods.hpp Conversions.hpp nucEOS.hpp
	$(CXX) $(CXXFLAGS) -c NumMethods.cpp Conversions.cpp nucEOS.cpp

quarkEOS.o: quarkEOS.hpp NumMethods.hpp
	$(CXX) $(CXXFLAGS) -c quarkEOS.cpp NumMethods.cpp

nucEOS.o: nucEOS.hpp quarkEOS.hpp NumMethods.hpp Conversions.hpp
	$(CXX) $(CXXFLAGS) -c nucEOS.cpp quarkEOS.cpp NumMethods.cpp Conversions.cpp

Conversions.o: Conversions.hpp
	$(CXX) $(CXXFLAGS) -c Conversions.cpp

MCMC.o: NumMethods.hpp nucEOS.hpp Conversions.hpp MCMC.hpp
	$(CXX) $(CXXFLAGS) -c NumMethods.cpp nucEOS.cpp Conversions.cpp MCMC.cpp

clean:
	rm -f main main.o NumMethods.o quarkEOS.o nucEOS.o Conversions.o MCMC.o