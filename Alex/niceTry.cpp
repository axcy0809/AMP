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


///////////////////////////////////////////////// TAS Lock

class LOCK 
{
    std::atomic<bool> state;
    
    public:
            // setup timer variables
        std::chrono::time_point<std::chrono::high_resolution_clock> start;
        std::chrono::time_point<std::chrono::high_resolution_clock> end;
        std::string file_name;
        double runtime = 0;
        double totruntime = 0;
        double HowFair = 0;
        double mean = 0;
        double var = 0;
        double tmp = 0;
        int tid;                                // thread ID
        long int counter = 0;                   // counter gets incremented in CS    

    LOCK()
    {
        state = false;
    }
    
    void TAS_lock(){
        while (state.exchange(true))
        {}
    }

    void TTAS_lock(){
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

    void run_TAS_lock(int numthreads, int iterations, int numofiter) 
    {
        // setup timer variables
        std::vector<int> Verteilung (numthreads,0);
        long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    
        for (int i = 0; i < numthreads; ++i)
            turns[std::max(i*8-1,0)] = 0;
        

        omp_set_num_threads(numthreads);        // setting number of threads
        for (int n = 0; n < numofiter; n++)
        {
            tmp = 0;
            start = std::chrono::high_resolution_clock::now();
            #pragma omp parallel private(tid) shared(counter)
            {		    
                tid = omp_get_thread_num();
                while(counter < iterations)
                {
                    TAS_lock();
                    try {
                        if(counter < iterations)
                        {
                            counter++;
                            turns[std::max(tid*8-1,0)]++;
                            std::cout << "Thread " << tid << " is in critical section" << std::endl;
                        }
                    }
                    catch (int j) {
                        //std::cout << "Some error occured while in CS" << std::endl;
                    }
                    unlock();
                }
                
            }
            end = std::chrono::high_resolution_clock::now();
            runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
            //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
            //std::cout << "counter " << counter << std::endl << std::endl;
            for (int i = 0; i < numthreads; ++i)
            {
                std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
                //std::cout << i <<  std::endl;
                Verteilung[i] = turns[std::max(i*8-1,0)];
            }

            mean = iterations/numthreads;

            for (int i = 0; i < numthreads; ++i)
            {
                tmp = tmp + (Verteilung[i] - mean) * (Verteilung[i] - mean);           
            }

            var = sqrt(tmp/numthreads);
            std::cout << var << std::endl;
            HowFair = HowFair + var;
            totruntime = totruntime + runtime;
            //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
            //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
            
        }
        file_name = "data/Max"+std::to_string(iterations)+"_"+std::to_string(numthreads)+"TAS_lock.csv";
            //return results;
        std::ofstream myfile;
        myfile.open (file_name);
        myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair" << "; \n";
        myfile << iterations << ";" << numthreads << ";" << totruntime/numofiter << ";" << HowFair/numofiter << ";";	  
        myfile << std::endl;

        myfile.close();
    }
    void run_TTAS_lock(int numthreads, int iterations, int numofiter) 
    {
        // setup timer variables                              // thread ID
        long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
        std::vector<int> Verteilung (numthreads,0);
        for (int i = 0; i < numthreads; ++i)
            turns[std::max(i*8-1,0)] = 0;
        

        omp_set_num_threads(numthreads);        // setting number of threads
        
        start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel private(tid) shared(counter)
        {		    
            tid = omp_get_thread_num();
            while(counter < iterations)
            {
                TTAS_lock();
                try 
                {
                    if(counter < iterations)
                    {
                        counter++;
                        turns[std::max(tid*8-1,0)]++;
                        std::cout << "Thread " << tid << " is in critical section" << std::endl;
                    }
                }
                catch (int j)
                {
                    std::cout << "Some error occured while in CS" << std::endl;
                }
                unlock();
            }
            
        }
        end = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        runtime = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        //std::cout << std::endl << "Program ran in TAS Lock with " << iterations << " iterations. Those are the results. " << std::endl << std::endl;
        //std::cout << "counter " << counter << std::endl << std::endl;
        for (int i = 0; i < numthreads; ++i)
        {
            std::cout << "turns[" << i << "] is " << turns[std::max(i*8-1,0)] << std::endl;
            //std::cout << i <<  std::endl;
            Verteilung[i] = turns[std::max(i*8-1,0)];
        }

        mean = iterations/numthreads;

        for (int i = 0; i < numthreads; ++i)
        {
            tmp = tmp + (Verteilung[i] - mean) * (Verteilung[i] - mean);           
        }

        var = sqrt(tmp/numthreads);
        std::cout << var << std::endl;
        HowFair = HowFair + var;
        totruntime = totruntime + runtime;
            //std::cout << std::endl << "runtime " << runtime << " s" << std::endl;
            //std::cout << std::endl << "iterations " << iterations << " numthreads " << results[iterations][0] << " runtime " << results[iterations][1] << " HowFair " << results[iterations][2] << std::endl;
            
        
        file_name = "data/Max"+std::to_string(iterations)+"_"+std::to_string(numthreads)+"TTAS_lock.csv";
                //return results;
        std::ofstream myfile;
        myfile.open (file_name);
        myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair" << "; \n";
        myfile << iterations << ";" << numthreads << ";" << totruntime/numofiter << ";" << HowFair/numofiter << ";";	  
        myfile << std::endl;

        myfile.close();

    }
};



///////////////////////////////////////////////////////////// main starts here //////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) 
{

    std::string file_name;
    std::vector<double> results;
    int mode = std::atoi(argv[1]);
    int numthreads = std::atoi(argv[2]);    // number of threads, command line input
    //long int maxiterations = std::atoi(argv[3]);              // number of iterations in CS
    long int iterations = std::atoi(argv[3]);
    results = std::vector<double>(3, 0);
    int numofiter = 50;
		
    LOCK Locker;
    //for (int iterations = 0; iterations < maxiterations; iterations = iterations + 1)
    //{
        if (mode == 1) 
        {
            Locker.run_TAS_lock(numthreads, iterations, numofiter);
        }
            
        if (mode == 2) 
        {
            Locker.run_TTAS_lock(numthreads, iterations, numofiter);
        }

    //}
    /*file_name = "data/"+std::to_string(numthreads)+"TAS_lock.csv";
    std::ofstream myfile;
    myfile.open (file_name);
    myfile << "NumberOfThreads" << ";" << "NumberOfIterations" << ";" << "RunTime" << ";" << "HowFair \n";
    myfile << numthreads << ";" << iterations << ";" << results[numthreads] << ";" << HowFair;	  
    myfile << std::endl;
    myfile.close();
    */
}

