
#include <iostream>
#include <iomanip>
#include <fstream>
#include <dirent.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <vector>
#include <chrono>
#include <map>
#include <sstream>

using namespace std;
int main(int argc, char *argv[])
{ 
	vector<int> test_1;
	test_1.push_back(5);
	test_1.push_back(3);

	vector<int> test_2;
	test_2 = test_1;

	test_2[1] = 10;

	cout << test_1[0] << " " << test_1[1] << endl;
	cout << test_2[0] << " " << test_2[1] << endl;
} 