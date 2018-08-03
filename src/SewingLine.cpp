#include "SewingLine.h"

namespace SeamStress
{
  void SewingLine<std::nullptr_t>::operator()(void *thread)
  {
    ulong w = *((ulong*)thread);
    pattern(w, nthreads, tapestry);
  }
  
  void SewingLine<std::nullptr_t>::sew(void (*design)(ulong, ulong, void*), void* Tapestry, long int n)
  {
    tapestry = Tapestry;
    pattern = design;
    nthreads = seamstresses->size();
    if(n > 0){nthreads = n;}
    quilt.resize(nthreads);
    bool own_thread = false;ulong which_own = 0;
    for(ulong i=0;i<nthreads;++i)
    {
      quilt[i] = i;
      if( (*seamstresses)[i] == nullptr ){own_thread = true;which_own = i;continue;}
      (*seamstresses)[i]->needle = this;
      (*seamstresses)[i]->thread = &(quilt[i]);
      (*seamstresses)[i]->sew();
    }
    if(own_thread == true){(*this)( &( quilt[which_own] ) );}
    for(unsigned long int i=0;i<nthreads;++i)
    {
      if( (*seamstresses)[i] == nullptr ){continue;}
      (*seamstresses)[i]->rest();
    }
  }
  
  void SewingLine<std::nullptr_t>::sewOpenly(void (*design)(ulong, ulong, void*), void* Tapestry, long int n)
  {
    tapestry = Tapestry;
    pattern = design;
    nthreads = seamstresses->size();
    if(n > 0){nthreads = n;}
    quilt.resize(nthreads);
    bool own_thread = false;ulong which_own = 0;
    for(ulong i=0;i<nthreads;++i)
    {
      quilt[i] = i;
      if( (*seamstresses)[i] == nullptr ){own_thread = true;which_own = i;continue;}
      (*seamstresses)[i]->needle = this;
      (*seamstresses)[i]->thread = &(quilt[i]);
      (*seamstresses)[i]->sew();
    }
    if(own_thread == true){(*this)( &( quilt[which_own] ) );}
  }
  
  void SewingLine<std::nullptr_t>::tieOff(long int n)
  {
    nthreads = seamstresses->size();
    if(n > 0){nthreads = n;}
    for(unsigned long int i=0;i<nthreads;++i)
    {
      if( (*seamstresses)[i] == nullptr ){continue;}
      (*seamstresses)[i]->rest();
    }
  }
}



