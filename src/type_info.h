#include <typeinfo>
#include <cxxabi.h>

#define DEMG(x) abi::__cxa_demangle(typeid(x).name(),0,0,&status)

#define VALUE_TEMPLATE_LSHIFT(x,y) std::cout << "Value of \"" << DEMG(x) << "\" " <<  x << " << " << y; std::cout << " is  " << (x << int_size<y>()) << "\n";
#define VALUE_TEMPLATE_RSHIFT(x,y) std::cout << "Value of \"" << DEMG(x) << "\" " <<  x << " >> " << y; std::cout << " is  " << (x >> int_size<y>()) << "\n";

#define VALUE_NAME_FUNC_INTARG1(a,F,b) F(a,b)
#define VALUE_NAME_SHIFT(x,y,OP) std::cout << "Value of " << DEMG(x) << " " << x << " "; std::cout << #OP << " " << y << " is  " << (x OP y) << "\n";
#define VALUE_NAME_FUNC(x, OP)	 std::cout << "Value of " << #OP << "(" << x << ") is " << OP(x) << "\n"; 


#define VALUE_NAME_OP(x,y,OP) \
	std::cout << "Value of " << DEMG(x) << " " << x << " " << #OP << " " << DEMG(y) << " "  << y << " is " << ( x OP y) << "\n"; 

#define VALUE_NAME_FUNC2(x,y, OP)					  \
	std::cout << "Value of " << #OP << "(" << DEMG(x) << x << "," << DEMG(y) << y << ") is " << OP(x,y) << "\n"; 


#define MIXED_VALUES(a,b,c)						\
  VALUE_NAME_OP( a, b, + );						\
  VALUE_NAME_OP( a, b, - );						\
  VALUE_NAME_OP( a, b, * );						\

#define THIS_MIXED_VALUES(a,b,c)						\
  VALUE_NAME_OP( a, b, += );						\
  VALUE_NAME_OP( a, b, -= );						\
  VALUE_NAME_OP( a, b, *= );						\


#define DMIXED_VALUES(a,b,c)						\
  VALUE_NAME_OP( a, b, + );						\
  VALUE_NAME_OP( b, a, + );						\
  VALUE_NAME_OP( a, b, - );						\
  VALUE_NAME_OP( b, a, - );						\
  VALUE_NAME_OP( b, a, * );						\
  VALUE_NAME_OP( a, b, * );						\

#define DMIXED_THIS_VALUES(a,b,c)						\
  VALUE_NAME_OP( a, b, += );						\
  VALUE_NAME_OP( b, a, += );						\
  VALUE_NAME_OP( a, b, -= );						\
  VALUE_NAME_OP( b, a, -= );						\
  VALUE_NAME_OP( b, a, *= );						\
  VALUE_NAME_OP( a, b, *= );						\


#define FUNC_VALUES(a,b,c)						\
  VALUE_NAME_FUNC( a, - );						\
  VALUE_NAME_SHIFT( a, b, >>); \
  VALUE_NAME_SHIFT( a, c, >>); \
  VALUE_NAME_SHIFT( a, b, <<); \
  VALUE_NAME_SHIFT( a, c, <<); 

#define FUNC2_VALUES(a,b,c) \
  VALUE_NAME_SHIFT( a, b, >>=); \
  VALUE_NAME_SHIFT( a, c, >>=); \
  VALUE_NAME_SHIFT( a, b, <<=); \
  VALUE_NAME_SHIFT( a, c, <<=); \

#define OTHER_FUNCS(a,b) \
  VALUE_NAME_FUNC2( a, b, pow); \
  VALUE_NAME_FUNC2( a, b, inverse);               \
  VALUE_NAME_FUNC( a, isPrime); 
