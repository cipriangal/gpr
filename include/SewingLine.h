#ifndef __SEWINGLINE__
#define __SEWINGLINE__

#include "Seamstress.h"
#include <tuple>
#include <functional>
#include <utility>


typedef unsigned long int ulong;

namespace SeamStress
{
  template <class TClass=std::nullptr_t>
  class SewingLine : public Needle
  {
    public:
      
      SewingLine<TClass>(TClass *glove, std::vector<Seamstress*> *ss) : seamstresses(ss), hand(glove)
      {
        
      }
      virtual ~SewingLine<TClass>(){}
      
      void operator()(void *thread)
      {
        ulong w = *((ulong*)thread);
        (hand->*pattern)(w, nthreads, tapestry);
      }
      
      std::vector<Seamstress*> *seamstresses;
      
      void sew(void (TClass::*design)(ulong, ulong, void*), void* Tapestry, long int n=-1)
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
      
      void sewOpenly(void (TClass::*design)(ulong, ulong, void*), void* Tapestry, long int n=-1)
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
      
      void tieOff(long int n=-1)
      {
        nthreads = seamstresses->size();
        if(n > 0){nthreads = n;}
        for(unsigned long int i=0;i<nthreads;++i)
        {
          if( (*seamstresses)[i] == nullptr ){continue;}
          (*seamstresses)[i]->rest();
        }
      }
      
    private:
      TClass *hand;
      void (TClass::*pattern)(ulong, ulong, void*);
      void* tapestry;
      ulong nthreads;
      std::vector<ulong> quilt;
  };
  
  template<>
  class SewingLine<std::nullptr_t> : public Needle
  {
    public:
      SewingLine<std::nullptr_t>(std::vector<Seamstress*> *ss) : seamstresses(ss) {}
      virtual ~SewingLine<std::nullptr_t>(){}
      std::vector<Seamstress*> *seamstresses;
      void operator()(void *thread);
      void sew(void (*design)(ulong, ulong, void*), void* Tapestry, long int n=-1);
      void sewOpenly(void (*design)(ulong, ulong, void*), void* Tapestry, long int n=-1);
      void tieOff(long int n=-1);
      
    private:
      template <typename T, typename... Ts>
      class Design
      {
      private:
        T f;
        std::tuple<Ts&...> args;
      public:
        Design(T F, Ts&... args) : f(F), args(std::tie(args...))
        {}
        
        template <typename... Args, std::size_t... Is>
        void func(ulong a, ulong b, std::tuple<Args...>& tup, std::integer_sequence<std::size_t, Is...>)
        {
          f(a, b, std::get<Is>(tup)...);
        }
        
        template <typename... Args>
        void func(ulong a, ulong b, std::tuple<Args...>& tup)
        {
          func(a, b, tup, std::make_index_sequence<sizeof...(Args)>{});
        }
        
        void act(ulong a, ulong b)
        {
          func(a, b, args);
        }
      };
      
    public:
      
      template <typename T, typename... ts>
      void Sew( T func, ts&... args )
      {
        Design<T, ts...> design( func, args... );
        tapestry = &design;
        pattern = []( ulong a, ulong b, void* A )
        {
          Design<T, ts...>* d = (Design<T, ts...>*)A;
          d->act(a, b);
        };
        nthreads = seamstresses->size();
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
      
    private:
      std::function< void(ulong, ulong, void*) > pattern;
      void* tapestry;
      ulong nthreads;
      std::vector<ulong> quilt;
  };
}

#endif

