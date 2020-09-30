all:
#	/usr/bin/env g++ -Wall -fPIC -o "simple-pendulum.o" "simple-pendulum.cpp"
	make docs
	make build-and-run
build-and-run:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.8" -I "/usr/include/python3.8/cpython" -I "/usr/lib/python3.8/site-packages/numpy/core/include" -lpython3.8 -o "RKF45.o" "RKF45.cpp"
	chmod +x RKF45.o
	./RKF45.o
build-docs:
	@doxygen ./Doxyfile
clean:
	rm RKF45.o
docs:
	make build-docs
	google-chrome-stable ./html/index.html