all:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -I ../headers -o "simplex.o" "simplex.cpp"
	chmod +x simplex.o
	./simplex.o

debug:
	/usr/bin/env g++ -g -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -I ../headers -o "simplex.o" "simplex.cpp"
	chmod +x simplex.o
	gdb simplex.o