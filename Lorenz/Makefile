# Requires support for C++ 11 (lorenz.cpp includes tuples)
# Clang 3.3 and later, as well as GCC 4.7 and later have at least experimental support for C++11
# 7m  52s is runtime with no -On flag
# 10m 25s is how long it takes to compile with -O2
# 4m  58s is how long it takes to compile with -O3
all:
	g++ -std=c++11 -O4 -pipe -Wall -fPIC -o "lorenz.o" "lorenz.cpp"
	chmod +x lorenz.o
	./lorenz.o
clean:
	rm lorenz.o
