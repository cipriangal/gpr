#ifndef __SEAMSTRESS__
#define __SEAMSTRESS__

#include <pthread.h>
#include "Needle.h"
#include <vector>
#include <atomic>

namespace SeamStress
{
  class Seamstress
  {
    public:
      Seamstress(bool aggro = false);
      Seamstress(const Seamstress &ss);
      ~Seamstress();
      
      static void *prepare(void *arg);
      
      void start();
      void stop();
      void sew();
      void rest();
      
      static void init_vector(unsigned long int N, std::vector<Seamstress> &vec);
      static std::vector<Seamstress*>* create_vector(unsigned long int N, unsigned long int naggro=0);
      
      
      Needle *needle;
      void *thread;
      pthread_mutex_t mutex;
      pthread_attr_t attr;
      
    private:
      std::atomic<bool> gotime;
      std::atomic<bool> end;
      std::atomic<bool> queue_end;
      std::atomic<bool> running;
      std::atomic<bool> started;
      pthread_t pthread;
      
      
      pthread_mutexattr_t mattr;
      pthread_cond_t cond;
      pthread_cond_t waitcond;
      
      bool aggressive;
  };
}


#endif
