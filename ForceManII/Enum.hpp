#pragma once

#include <string>
#include <set>


/** \brief A class to implement intelligent enumerations, based on strings
 *
 *  I stole this class from drdobbs.com.  I've renamed stuff and
 *  made a few modifications so that it reflects my coding style, but
 *  really all credit goes to that website for the pattern and implementation.
 *
 *  Anyways, this class is geared at making what I will call intelligent
 *  enumerations.  Basically, these enumerations support unbreakable iteration
 *  (probably not a big concern for most users) which means even if
 *  the enums are not sequential, e.g. Enum1=1, Enum2=3,etc. you can still
 *  iterate over them.  More important for my liking, they have the desirable
 *  property of being scoped like the C++11 strong enums, and they know how to
 *  print their values!!!!
 *
 *  Note that these are fundamentally strings so there is some overhead
 *  associated with the comparisons that occur.
 *
 *  To declare a set of enums called MyEnums, with values: Value1,Value2, and
 *  Value3 do in a header file MyEnums.hpp:
 *  \code
 *  class MyEnums
 *       : public Enumeration<MyEnums>{
 *     private:
 *        MyEnums(const std::string& Name):Enumeration<MyEnums>(Name){}
 *     public:
 *        static const MyEnums Value1;
 *        static const MyEnums Value2;
 *        static const MyEnums Value3;
 *  };
 *  \endcode
 *
 *  Then in a source file, MyEnums.cpp:
 *  \code
 *  #include "MyEnums.hpp"
 *  std::set<const Enumeration<MyEnums>*> Enumeration<MyEnums>::Enums_;
 *  const MyEnums MyEnums::Value1("String to print for Value1's value");
 *  const MyEnums MyEnums::Value2("String to print for Value2's value");
 *  const MyEnums MyEnums::Value3("String to print for Value3's value");
 *  \endcode
 *
 *  The first line instantiates an std::set to hold all of the enums
 *  where as the rest of the lines instantiate the actual enums.
 *
 *  Then in code where you want to use the enums:
 *  \code
 *  #include "MyEnums.hpp"
 *
 *  void FunctionThatTakesEnum(const MyEnums& AnEnum){
 *     if(AnEnum==MyEnums::Value1)//Do stuff for Value1
 *     else if(AnEnum==MyEnums::Value2)//Do stuff for Value2
 *     else //Do stuff for Value3
 *  }
 *  \endcode
 */
template <class T>
class Enumeration{
   private:
      typedef Enumeration<T> MyType_t;
      typedef std::set<const MyType_t*> Set_t;

      ///The value this enumeration holds
      std::string Value_;

      ///A list of all enumerations, common to all Enumeration<T>
      static Set_t Enums_;

   protected:
      /// Constructor, protected cause you won't ever call it outside derived class
      explicit Enumeration(const std::string& InValue):Value_(InValue){
         Enums_.insert(this);
      }

   public:
      bool operator<(const MyType_t& Other)const{return Value()<Other.Value();}
      bool operator>(const MyType_t& Other)const{return Other<(*this);}
      bool operator>=(const MyType_t& Other)const{return !operator<(Other);}
      bool operator<=(const MyType_t& Other)const{return !operator>(Other);}
      bool operator==(const MyType_t& Other)const{
         return (!operator<(Other)&&!operator>(Other));
      }
      bool operator==(const std::string& other)const{
         if(IsValid(other))
            return (*this)==MyType_t(other);
         return false;
      }
      bool operator!=(const MyType_t& Other)const{return !((*this)==Other);}
      ///Returns the value held by this Enumeration
      const std::string& Value() const{return Value_;}

      operator std::string()const{return Value();}

      static bool IsValid(const std::string& Value){return
            Enums_.count(MyType_t(Value))==1; }

      // Number of elements
      static size_t size(){return Enums_.size();}

      //@{
      ///Iterating over enumerations
      ///Returns the minimum value, for iterating in a for loop w/o iterators
      const std::string& Min()const{return (*Enums_.begin())->Value_;}
      ///Returns the last enumeration value, for iterating w/o iterators
      const std::string& Max()const{return (*Enums_.rbegin())->Value_;}
      typedef typename Set_t::const_iterator const_iterator;
      static const_iterator begin(){ return Enums_.begin(); }
      static const_iterator end(){ return Enums_.end(); }
      //@}
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os,const Enumeration<T>& Enum){
   return os<<(std::string)(Enum);
}

template<typename T>
inline bool operator==(const std::string& lhs,const Enumeration<T>& rhs){
   return rhs==lhs;
}

//template<typename T>
//std::set<const Enumeration<T>* > Enumeration<T>::Enums_;

}}//end namespaces
