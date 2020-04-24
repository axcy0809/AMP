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

///////////////////////////////////////////////// TAS Lock

class TAS_lock {
    std::atomic<bool> state;
    
    public:

    TAS_lock(){
        state = false;
    }
    
    void lock(){
        while (state.exchange(true))
        {}
    }

    void unlock(){
        state.exchange(false);
    }

};
        
//////////////////////////////////////////////// TTAS Lock

class TTAS_lock {
    std::atomic<bool> state;
    
    public:

    TTAS_lock(){
        state = false;
    }
    
    void lock(){
        while (true) {
            while (state.load()) 
            {}
            if (!state.exchange(true))
                return;
        }
    }

    void unlock(){
        state.exchange(false);
    }

};

/////////////////////////////////////////////// TicketLock

class Ticket_lock {
    std::atomic<int> ticket;
    volatile int served;

    public:
    
    Ticket_lock() {
        ticket = 0;
        served = 0;
    }

    void lock() {
        int next = ticket.fetch_add(1);
        while (served < next)
        {}
    }

    void unlock() {
        served++;
    }

};


///////////////////////////////////////////////////////////// main starts here //////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) 
{
    // setup timer variables
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    double runtime;

    int tid;                                // thread ID
    int numthreads = std::atoi(argv[1]);    // number of threads, command line input
    long int iterations = 1e2;              // number of iterations in CS
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    Ticket_lock mylock;                         // instantiate the lock we will use

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		
        tid = omp_get_thread_num();

        while(counter < iterations)
        {
            mylock.lock();
            // Critical Section
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            mylock.unlock();
        }
    }
    end = std::chrono::high_resolution_clock::now();
    runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << "Program ran with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
    std::cout << "counter " << counter << std::endl << std::endl;
    for (int i = 0; i < numthreads; ++i)
        std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
    std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
				
}


