doublePendulum:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "doublePendulum.o" "doublePendulum.cpp"
	chmod +x doublePendulum.o
	./doublePendulum.o

Lorenz:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "Lorenz.o" "Lorenz.cpp"
	chmod +x Lorenz.o
	./Lorenz.o

simplePendulum:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "simplePendulum.o" "simplePendulum.cpp"
	chmod +x simplePendulum.o
	./simplePendulum.o

SIR:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "SIR.o" "SIR.cpp"
	chmod +x SIR.o
	./SIR.o

debug-doublePendulum:
	/usr/bin/env g++ -g -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "doublePendulum.o" "doublePendulum.cpp"
	chmod +x doublePendulum.o
	gdb doublePendulum.o

debug-Lorenz:
	/usr/bin/env g++ -g -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "Lorenz.o" "Lorenz.cpp"
	chmod +x Lorenz.o
	gdb Lorenz.o

debug-simplePendulum:
	/usr/bin/env g++ -g -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "simplePendulum.o" "simplePendulum.cpp"
	chmod +x simplePendulum.o
	gdb simplePendulum.o

debug-SIR:
	/usr/bin/env g++ -g -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "SIR.o" "SIR.cpp"
	chmod +x SIR.o
	gdb SIR.o
