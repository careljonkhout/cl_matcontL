#define SUNDIALS_DOUBLE_PRECISION 1
//#define SUNDIALS_EXTENDED_PRECISION 1
#define HIGHER_PRECISION __float128
// comment line 1 and uncomment line 2 to enable higher precision
// (i.e. an 80 or 128 bit double). Extended precision uses more memory ( 
// storing an extended precision double probably costs uses 16 bytes ), and more
// cpu time
//
// valid values for HIGHER_PRECISION include
//         "long double"
//         "__float128" ( for gcc )
//
// Note that the implementation of long double is compiler specific. 