#include <.vecOpsFunctional.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
using namespace std::chrono;
std::vector<double> vec1;
vec1 = {1.5, 2.3, 4.5, 6.7, 9.2, 2.23, 5.5};
std::vector<double> vec2;
vec2 = {1.5, 2.3, 4.5, 6.7, 9.2, 2.23, 8.5};
std::vector<double> vec3;
vec3 = {1.5, 2.3, 4.5, 3.7, 6.2, 2.23, 8.5};

auto start1 = high_resolution_clock::now();
vecAdd({vec1, vec2});
auto end1 = high_resolution_clock::now();

auto start2 = high_resolution_clock::now();
vecAdd2(vec1, vec2);
auto end2 = high_resolution_clock::now();

auto start3 = high_resolution_clock::now();
vecAdd3(vec1, vec2, vec3);
auto end3 = high_resolution_clock::now();

auto start4 = high_resolution_clock::now();
vecAdd({vec1, vec2, vec3});
auto end4 = high_resolution_clock::now();

auto startTrans12 = high_resolution_clock::now();
vecAddTrans({vec1, vec2});
auto endTrans12 = high_resolution_clock::now();

auto startTrans123 = high_resolution_clock::now();
vecAddTrans({vec1, vec2, vec3});
auto endTrans123 = high_resolution_clock::now();

auto duration1 = duration_cast<milliseconds>(end1 - start1);
std::cout << "vecAdd({vec1, vec2}) took " << duration1.count() << " ms to run";

auto duration2 = duration_cast<milliseconds>(end2 - start2);
std::cout << "vecAdd2(vec1, vec2) took " << duration2.count() << " ms to run";

auto duration4 = duration_cast<milliseconds>(end4 - start4);
std::cout << "vecAdd({vec1, vec2, vec3}) took " << duration4.count() << " ms to run";

auto duration3 = duration_cast<milliseconds>(end3 - start3);
std::cout << "vecAdd3(vec1, vec2, vec3) took " << duration3.count() << " ms to run";

auto duration12 = duration_cast<milliseconds>(endTrans12 - startTrans12);
std::cout << "vecAddTrans({vec1, vec2}) took " << duration12.count() << " ms to run";

auto duration123 = duration_cast<milliseconds>(endTrans123 - startTrans123);
std::cout << "vecAddTrans({vec1, vec2, vec3}) took " << duration123.count() << " ms to run";