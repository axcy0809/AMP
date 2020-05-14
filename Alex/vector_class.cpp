#include <assert.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <omp.h>
#include <atomic>
#include <chrono>
#include <stdlib.h>
#include <thread>
#include <string>
#include <fstream>

using namespace std; 
class Geeks 
{ 
    // Access specifier 
    public: 
  
    // Data Members 
    string geekname; 
  
    // Member Functions() 
    void printname() 
    { 
       cout << "Geekname is: " << geekname; 
    } 
}; 
  
int main() { 
  
    // Declare an object of class geeks 
    Geeks obj1; 
  
    // accessing data member 
    obj1.geekname = "Abhi"; 
  
    // accessing member function 
    obj1.printname(); 
    return 0; 
} 
