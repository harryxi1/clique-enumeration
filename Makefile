 output: colour.o pivoter.o utils.o
		g++ -std=c++0x -Wall colour.o pivoter.o utils.o -o output

 colour.o: colour.cpp
		g++ -c colour.cpp

 pivoter.o: pivoter.cpp
		g++ -c pivoter.cpp

 utils.o: utils.cpp
		g++ -c utils.cpp

 clean:
		rm *.o output