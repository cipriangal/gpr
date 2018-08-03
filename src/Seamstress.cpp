#include "Seamstress.h"
#include <iostream>
#include <unistd.h>
using namespace std;


namespace SeamStress
{
  Seamstress::Seamstress( bool aggro ) : aggressive(aggro)
  {
    gotime = false;
    end = false;
    queue_end = false;
    running = false;
    started = false;
    pthread_mutexattr_init(&mattr);
    pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_NORMAL);
    pthread_mutex_init(&mutex, &mattr);
    pthread_cond_init(&cond, NULL);
    pthread_cond_init(&waitcond, NULL);
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  }


  Seamstress::Seamstress(const Seamstress &ss)
  {
    gotime = false;
    end = false;
    queue_end = false;
    running = false;
    started = false;
    pthread_mutexattr_init(&mattr);
    pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_NORMAL);
    pthread_mutex_init(&mutex, &mattr);
    pthread_cond_init(&cond, NULL);
    pthread_cond_init(&waitcond, NULL);
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  }


  Seamstress::~Seamstress()
  {
    pthread_mutex_lock(&mutex);
    while(started==true)
    {
      pthread_cond_wait(&waitcond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    
    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&mutex);
    pthread_mutexattr_destroy(&mattr);
    pthread_cond_destroy(&cond);
    pthread_cond_destroy(&waitcond);
  }


  void *Seamstress::prepare(void *arg)
  {
    Seamstress *seamstress = (Seamstress*)arg;
    
    if(seamstress->aggressive == false)
    {
      while(true)
      {
        pthread_mutex_lock(&(seamstress->mutex));
        if((seamstress->end)==true){pthread_mutex_unlock(&(seamstress->mutex));break;}
        if(seamstress->queue_end == true)
        {
          seamstress->end = true;
          pthread_mutex_unlock(&(seamstress->mutex));
          continue;
        }
        else
        {
          if((seamstress->end)==true)
          {
            pthread_mutex_unlock(&(seamstress->mutex));
            break;
          }
          while(seamstress->gotime==false)
          {
            pthread_cond_wait(&(seamstress->cond), &(seamstress->mutex));
          }
          if(seamstress->queue_end == true)
          {
            seamstress->end = true;
            pthread_mutex_unlock(&(seamstress->mutex));
            continue;
          }
        }
        
        seamstress->gotime = false;
        (*(seamstress->needle))(seamstress->thread);
        
        seamstress->running = false; 
        if(seamstress->queue_end==true){seamstress->end=true;}
        pthread_mutex_unlock(&(seamstress->mutex));
        pthread_cond_signal(&(seamstress->waitcond));
      }
      pthread_mutex_lock(&(seamstress->mutex));
      seamstress->started = false;
      seamstress->running = false;
      pthread_mutex_unlock(&(seamstress->mutex));
      pthread_cond_signal(&(seamstress->waitcond));
      pthread_exit(NULL);
    }
    else
    {
      while(true)
      {
        while( true )
        {
          if( (seamstress->queue_end).load(memory_order_relaxed) == true ){break;}
          if( (seamstress->gotime).load(memory_order_relaxed) == true ){break;}
        }
        
        if( (seamstress->queue_end).load(memory_order_relaxed) == true ){break;}
        
        seamstress->gotime = false;
        (*(seamstress->needle))(seamstress->thread);
        
        seamstress->running = false;
        
      }
      
      pthread_mutex_lock(&(seamstress->mutex));
      seamstress->started = false;
      seamstress->running = false;
      pthread_mutex_unlock(&(seamstress->mutex));
      pthread_cond_signal(&(seamstress->waitcond));
      
      
      pthread_exit(NULL);
    }
  }


  void Seamstress::start()
  {
    pthread_mutex_lock(&mutex);
    if(started==false)
    {
      end = false;
      gotime = false;
      queue_end = false;
      running = false;
      started = true;
      pthread_mutex_unlock(&mutex);
      pthread_create(&pthread, &attr, &Seamstress::prepare, this);
    }
    else
    {
      pthread_mutex_unlock(&mutex);
    }
  }


  void Seamstress::stop()
  {
    if(aggressive == false)
    {
      pthread_mutex_lock(&mutex);
      queue_end = true;
      gotime = true;
      pthread_mutex_unlock(&mutex);
      pthread_cond_signal(&cond);
    }
    else
    {
      queue_end = true;
    }
  }


  void Seamstress::sew()
  {
    if(aggressive == false)
    {
      pthread_mutex_lock(&mutex);
      gotime = true;
      running = true;
      pthread_mutex_unlock(&mutex);
      pthread_cond_signal(&cond);
    }
    else
    {
      running = true;
      gotime = true;
    }
  }


  void Seamstress::rest()
  {
    if(aggressive == false)
    {
      pthread_mutex_lock(&mutex);
      while(running==true)
      {
        pthread_cond_wait(&waitcond, &mutex);
      }
      pthread_mutex_unlock(&mutex);
    }
    else
    {
      while( true )
      {
        if(running == false){break;}
      }
    }
  }


  void Seamstress::init_vector(unsigned long int N, vector<Seamstress> &vec)
  {
    vec.clear();
    Seamstress st;
    
    for(unsigned int i=0;i<N;i++)
    {
      vec.push_back(st);
    }
    for(unsigned int i=0;i<N;i++)
    {
      vec[i].start();
    }
  }
  
  
  vector<Seamstress*>* Seamstress::create_vector(unsigned long int N, unsigned long int naggro)
  {
    vector<Seamstress*>* vec = new vector<Seamstress*>();
    
//     for(unsigned int i=0;i<naggro;i++){vec->push_back(new Seamstress(true));}
//     for(unsigned int i=naggro;i<N;i++){vec->push_back(new Seamstress(false));}
    
    for(unsigned int i=0;i<(N-naggro);i++){vec->push_back(new Seamstress(false));}
    for(unsigned int i=(N-naggro);i<N;i++){vec->push_back(new Seamstress(true));}
    
    
    for(unsigned long int i=0;i<N;i++){(*vec)[i]->start();}
    return vec;
  }
}


