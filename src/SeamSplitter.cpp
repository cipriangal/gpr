#include "SeamSplitter.h"

namespace SeamStress
{
  void SeamSplitter<std::nullptr_t>::operator()(void *thread)
  {
    split_struct* split = (split_struct*)thread;
    (split->pattern)(&half_1_ss, split->tap);
  }
  
  void SeamSplitter<std::nullptr_t>::split( void (*design1)(std::vector<Seamstress*>*, void*), void* Tapestry1, void (*design2)(std::vector<Seamstress*>*, void*), void* Tapestry2 )
  {
    if( seamstresses->size() > 1 )
    {
      headseamstress = (*seamstresses)[0];
      half_1_ss.clear();
      for(ulong i=1, loopsizei=((seamstresses->size())/2);i<loopsizei;i+=1){half_1_ss.push_back( (*seamstresses)[i] );}
      half_1_ss.push_back(nullptr);
      split1.tap = Tapestry1;split1.pattern = design1;
      headseamstress->needle = this;headseamstress->thread = (void*)(&split1);
      headseamstress->sew();
      
      half_2_ss.clear();
      for(ulong i=((seamstresses->size())/2), loopsizei=(seamstresses->size());i<loopsizei;i+=1){half_2_ss.push_back( (*seamstresses)[i] );}
      (*design2)( &half_2_ss, Tapestry2 );
      headseamstress->rest();
    }
    else
    {
      (*design1)( seamstresses, Tapestry1 );
      (*design2)( seamstresses, Tapestry2 );
    }
  }
}
