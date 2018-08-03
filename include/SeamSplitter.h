#ifndef __SEAMSPLITTER__
#define __SEAMSPLITTER__

#include "Seamstress.h"
#include <tuple>
#include <functional>
#include <utility>

typedef unsigned long int ulong;

namespace SeamStress
{
  template <class TClass=std::nullptr_t>
  class SeamSplitter : public Needle
  {
    public:
      std::vector<Seamstress*> *seamstresses;
      SeamSplitter<TClass>(TClass* glove, std::vector<Seamstress*> *ss) : seamstresses(ss), hand(glove) {}
      virtual ~SeamSplitter<TClass>(){}
      
      void operator()(void *thread)
      {
        split_struct* split = (split_struct*)thread;
        (hand->*(split->pattern))(&half_1_ss, split->tap);
      }
      
      void split( void (TClass::*design1)(std::vector<Seamstress*>*, void*), void* Tapestry1, void (TClass::*design2)(std::vector<Seamstress*>*, void*), void* Tapestry2 )
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
          (hand->*design2)( &half_2_ss, Tapestry2 );
          headseamstress->rest();
        }
        else
        {
          (hand->*design1)( seamstresses, Tapestry1 );
          (hand->*design2)( seamstresses, Tapestry2 );
        }
      }
      
    private:
      
      struct split_struct
      {
        void (TClass::*pattern)(std::vector<Seamstress*>*, void*);
        void* tap;
      };
      TClass* hand;
      Seamstress* headseamstress;
      std::vector<Seamstress*> half_1_ss;
      std::vector<Seamstress*> half_2_ss;
      split_struct split1;
  };
  
  template<>
  class SeamSplitter<std::nullptr_t> : public Needle
  {
    public:
      SeamSplitter<std::nullptr_t>(std::vector<Seamstress*> *ss) : seamstresses(ss) {}
      virtual ~SeamSplitter<std::nullptr_t>(){}
      std::vector<Seamstress*> *seamstresses;
      void operator()(void *thread);
      void split( void (*design1)(std::vector<Seamstress*>*, void*), void* Tapestry1, void (*design2)(std::vector<Seamstress*>*, void*), void* Tapestry2 );
      
      template <typename... Ts1, typename... Ts2>
      void Split( std::tuple<Ts1&...>&& t1, std::tuple<Ts2&...>&& t2 )
      {
        split2( t1, shift_sequence(std::make_index_sequence< (sizeof...(Ts1))-1 >{}), t2, shift_sequence(std::make_index_sequence< (sizeof...(Ts2))-1 >{}) );
      }
      
    private:
      template <typename T, typename... Ts>
      class Design
      {
      private:
        T f;
        std::tuple<Ts&...> args;
      public:
        Design(T F, Ts&... Args) : f(F), args( tie(Args...) ) {}
        
        template <typename... Args, size_t... Is>
        void func(std::vector<Seamstress*>* a, std::tuple<Args...>& tup, std::integer_sequence<std::size_t, Is...>)
        {
          f(a, std::get<Is>(tup)...);
        }
        
        template <typename... Args>
        void func(std::vector<Seamstress*>* a, std::tuple<Args...>& tup)
        {
          func(a, tup, std::make_index_sequence<sizeof...(Args)>{});
        }
        
        void act(std::vector<Seamstress*>* a)
        {
          func(a, args);
        }
      };
      
      template <typename T, typename... Ts>
      Design<T, Ts...> make_design( T F, Ts&... Args )
      {
        return Design<T, Ts...>(F, Args...);
      }
      
      template <typename... Ts1, typename... Ts2>
      void split3( Design<Ts1...>&& d1, Design<Ts2...>&& d2 )
      {
        if( seamstresses->size() > 1 )
        {
          split1.tap = (void*)(&d1);
          split1.pattern = [](std::vector<Seamstress*>* s, void* A)
          {
            Design<Ts1...>* d = (Design<Ts1...>*)A;
            d->act(s);
          };
          headseamstress = (*seamstresses)[0];
          half_1_ss.clear();
          for(ulong i=1, loopsizei=((seamstresses->size())/2);i<loopsizei;i+=1){half_1_ss.push_back( (*seamstresses)[i] );}
          half_1_ss.push_back(nullptr);
          headseamstress->needle = this;headseamstress->thread = (void*)(&split1);
          headseamstress->sew();
          
          half_2_ss.clear();
          for(ulong i=((seamstresses->size())/2), loopsizei=(seamstresses->size());i<loopsizei;i+=1){half_2_ss.push_back( (*seamstresses)[i] );}
          d2.act( &half_2_ss );
          headseamstress->rest();
        }
        else
        {
          d1.act( seamstresses );
          d2.act( seamstresses );
        }
      }
      
      template <std::size_t... Is>
      std::integer_sequence<std::size_t, 1 + Is...> shift_sequence( std::integer_sequence<std::size_t, Is...> )
      {
        return std::integer_sequence<std::size_t, 1 + Is...>{};
      }
      
      template <typename... Ts1, std::size_t... Is1, typename... Ts2, std::size_t... Is2>
      void split2( std::tuple<Ts1&...>& t1, std::integer_sequence<std::size_t, Is1...>, std::tuple<Ts2&...>& t2, std::integer_sequence<std::size_t, Is2...> )
      {
        split3( make_design( std::get<0>(t1),  (std::get<Is1>(t1))... ), make_design( std::get<0>(t2),  (std::get<Is2>(t2))... ) );
      }
      
      
      struct split_struct
      {
        std::function< void(std::vector<Seamstress*>*, void*) > pattern;
        void* tap;
      };
      Seamstress* headseamstress;
      std::vector<Seamstress*> half_1_ss;
      std::vector<Seamstress*> half_2_ss;
      split_struct split1;
  };
}


#endif
