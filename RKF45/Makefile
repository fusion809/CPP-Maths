simplePendulum:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "simplePendulum.o" "simplePendulum.cpp"
	chmod +x simplePendulum.o
	./simplePendulum.o

Lorenz:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "Lorenz.o" "Lorenz.cpp"
	chmod +x Lorenz.o
	./Lorenz.o

debug:
	/usr/bin/env g++ -g -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -I . -o "RKF45.o" "RKF45.cpp"
	chmod +x RKF45.o
	gdb RKF45.o
