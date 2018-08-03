#ifndef __PINCUSHION__
#define __PINCUSHION__

#include "Seamstress.h"


namespace SeamStress
{
  template <class TClass>
  class Pincushion : public Needle
  {
    public:
      Pincushion<TClass>(TClass *glove, std::vector<Seamstress*> *ss)
      {
        hand = glove;
        seamstresses = ss;
      }
      virtual ~Pincushion<TClass>(){}
      
      
      void operator()(void *thread)
      {
        (hand->*pattern)(thread);
      }
      
      
      void sew(void (TClass::*design)(void*), std::vector<void*> threads)
      {
        if(threads.size()==0){return;}
        pattern = design;
        if(seamstresses==NULL){(*this)(threads[0]);return;}
        unsigned long int ntot = threads.size();
        if(seamstresses->size() < ntot){ntot = seamstresses->size();}
        for(unsigned long int i=0;i<ntot;i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(threads[i]);continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = threads[i];
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<ntot;i++)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewOpenly(void (TClass::*design)(void*), std::vector<void*> threads)
      {
        if(threads.size()==0){return;}
        pattern = design;
        unsigned long int ntot = threads.size();
        if(seamstresses->size() < ntot){ntot = seamstresses->size();}
        for(unsigned long int i=0;i<ntot;i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(threads[i]);continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = threads[i];
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void sewSoftly(void (TClass::*design)(void*))
      {
        pattern = design;
        if(seamstresses==NULL){(*this)(NULL);return;}
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(NULL);continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewSoftly(void (TClass::*design)(void*), unsigned long int N)
      {
        pattern = design;
        if(seamstresses==NULL){(*this)(NULL);return;}
        for(unsigned long int i=0;i<N;i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(NULL);continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<N;i++)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewStraight(void (TClass::*design)(void*))
      {
        pattern = design;
        if(seamstresses==NULL){unsigned long int val=0;(*this)(&val);return;}
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<seamstresses->size();i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(&(quilt[i]));continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewStraight(void (TClass::*design)(void*), unsigned long int N)
      {
        if(N==0){return;}
        pattern = design;
        if(seamstresses==NULL){unsigned long int val=0;(*this)(&val);return;}
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<N;i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<N;i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(&(quilt[i]));continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<N;i++)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewOpenlyStraight(void (TClass::*design)(void*))
      {
        pattern = design;
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<seamstresses->size();i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(&(quilt[i]));continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void sewOpenlyStraight(void (TClass::*design)(void*), unsigned long int N)
      {
        pattern = design;
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<N;i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<N;i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(&(quilt[i]));continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void sewOpenlySoft(void (TClass::*design)(void*))
      {
        pattern = design;
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(NULL);continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void sewOpenlySoft(void (TClass::*design)(void*), unsigned long int N)
      {
        pattern = design;
        for(unsigned long int i=0;i<N;i++)
        {
          if( (*seamstresses)[i] == nullptr ){(hand->*pattern)(NULL);continue;}
          pthread_mutex_lock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          pthread_mutex_unlock(&((*seamstresses)[i]->mutex));
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void tieOff()
      {
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void tieOff(unsigned long int N)
      {
        for(unsigned long int i=0;i<N;i++)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
      
      std::vector<Seamstress*> *seamstresses;
      TClass *hand;
      void (TClass::*pattern)(void*);
  };
}


#endif
