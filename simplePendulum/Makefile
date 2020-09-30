all:
#	/usr/bin/env g++ -Wall -fPIC -o "simple-pendulum.o" "simple-pendulum.cpp"
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -o "RKF45.o" "RKF45.cpp"
	chmod +x RKF45.o
	./RKF45.o
clean:
	rm RKF45.o
docs:
	@doxygen ./Doxyfile
	google-chrome-stable ./html/index.html