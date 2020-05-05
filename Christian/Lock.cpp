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

/////////////////////////////////////////////// Ticket Lock

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

///////////////////////////////////////////////// Array Lock
/*
class Array_lock {
    bool* flag;
    std::atomic<int> tail;
    std::atomic<int>* mySlot;
    int numthreads;

    public:

    Array_lock(int n) : flag(new bool[n]), mySlot(new std::atomic<int>[n]) {    ///////////////////////////// false sharing noch korrigieren
        for (int i = 0; i < n; ++i)
            flag[i] = false;
        flag[0] = true;
        tail = 0;
        numthreads = n;
    }

    void lock(int tid) {
        mySlot[tid] = tail.fetch_add(1)%numthreads;
        while (!flag[mySlot[tid]])
        {}
    }

    void unlock(int tid) {
        flag[mySlot[tid]] = false;
        flag[(mySlot[tid]+1)%numthreads] = true;
    }

};*/

class Array_lock {
    std::atomic<bool>* flag;
    std::atomic<int> tail;
    int numthreads;

    public:

    Array_lock(int n) : flag(new std::atomic<bool>[n]) {    
        for (int i = 0; i < n; ++i)
            flag[i] = false;
        flag[0] = true;
        tail = 0;
        numthreads = n;
    }

    void lock(int* mySlot) {
        *mySlot = tail.fetch_add(1)%numthreads;
        while (!flag[*mySlot])
        {}
    }

    void unlock(int* mySlot) {
        flag[*mySlot] = false;
        flag[(*mySlot+1)%numthreads] = true;
    }

};


///////////////////////////////////////////////////////////// CLH Lock

class QNode {

    public:
    std::atomic<bool> locked;
    QNode* pred;
    
    QNode() {
        locked = false;
        pred = nullptr;
    }   
};
    
class CLH_lock {

    std::atomic<QNode*> tail;
    
    public:
    
    CLH_lock() {
        tail = new QNode;
    }

    void lock(QNode** pointerToNode) {
        QNode* node = new QNode;
        *pointerToNode = node;
        node->locked = true;
        node->pred = std::atomic_exchange(&tail,node);
        while (node->pred->locked) 
        {}
    }

    void unlock(QNode* node) {  
        delete node->pred;
        node->locked = false;
    }

};

///////////////////////////////////////////////////////////////////////// MCS Lock

class Node {
    
    public:
        std::atomic<bool> locked;
        Node* pred;
        Node* next;

    Node() {
        pred = nullptr;
        next = nullptr;
        locked = false;
    }
};

class MCS_lock {
    
    std::atomic<Node*> tail;

    public:

    MCS_lock() {
        tail = nullptr;
    }

    void lock(Node* node, int tid) {
        Node* my = node;
        Node* pred = std::atomic_exchange(&tail,my);
        if (pred != nullptr) {
            my->locked = true;
            pred->next = my;
            while (my->locked)
            {}
        }
    }

    void unlock(Node* node, int tid) {
        Node* my = node;
        if (my->next == nullptr) {
            Node* compare = tail.exchange(nullptr);
            if (compare == my) {
                return;
            }
            my = node;
            while (my->next == nullptr)
            {}
        }
        my->next->locked = false;
        my->next =  nullptr;
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
    long int iterations = std::atoi(argv[2]);              // number of iterations in CS
    long int counter = 0;                   // counter gets incremented in CS
    long int turns[(numthreads-1)*8+1];     // keeps count of how often a thread got the CS, long has 8 bytes, so write one value every 8 slots to avoid false sharing
    //TTAS_lock mylock;
    //Array_lock mylock(numthreads);          // instantiate the lock we will use
    //CLH_lock mylock;
    MCS_lock mylock;

    for (int i = 0; i < numthreads; ++i)
        turns[std::max(i*8-1,0)] = 0;
    

    omp_set_num_threads(numthreads);        // setting number of threads
	
    start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel private(tid) shared(counter)
	{		
        //thread_local int mySlot;
        //thread_local QNode* pointerToNode;
        thread_local Node* my = new Node();
        tid = omp_get_thread_num();

        while(counter < iterations)
        {
            //mylock.lock();
            //mylock.lock(&pointerToNode);
            mylock.lock(my, tid);
            //mylock.lock(&mySlot);
            // Critical Section
            try {
                if(counter < iterations)
                {
                    counter++;
                    turns[std::max(tid*8-1,0)]++;
                    std::cout << "Thread " << tid << " is in critical section" << std::endl;
                }
            }
            catch (int j) {
                std::cout << "Some error occured while in CS" << std::endl;
            }
            //mylock.unlock();
            //mylock.unlock(pointerToNode);
            mylock.unlock(my, tid);
            //mylock.unlock(&mySlot);
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


