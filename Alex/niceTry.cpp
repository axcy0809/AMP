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
class Locks
{

};

class TAS : public Locks
{
    std::atomic<bool> state;
    
    public:

    TAS()
    {
        state = false;
    }
    
    void lock()
    {
        while (state.exchange(true))
        {}
    }

    void unlock()
    {
        state.exchange(false);
    }

};