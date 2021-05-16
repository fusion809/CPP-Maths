all:
#	/usr/bin/env g++ -Wall -fPIC -o "simple-pendulum.o" "simple-pendulum.cpp"
	make docs
	make build-and-run
build:
	/usr/bin/env g++ -Wall -fPIC -I "/usr/include/python3.9" -I "/usr/include/python3.9/cpython" -I "/usr/lib/python3.9/site-packages/numpy/core/include" -lpython3.9 -o "RKF45.o" "RKF45.cpp"
build-and-run:
	make build
	chmod +x RKF45.o
	./RKF45.o
build-docs:
	@doxygen ./Doxyfile
clean:
	rm RKF45.o
docs:
	make build-docs
	google-chrome-stable ./html/index.html