#define EXTENDED_PREC_EPSILON 1.084202172485504e-19
#define DOUBLE_DOUBLE_EPSILON 1.232595164407831e-32
//#define QUAD_EPSILON          1.925929944387236e-34


#define SUNDIALS_DOUBLE_PRECISION 1
//#define SUNDIALS_EXTENDED_PRECISION 1
#define HIGHER_PRECISION __float128
#define HIGHER_PRECISION_EPSILON FLT128_EPSILON


// comment line 6 and uncomment line 7 and 8 to enable higher precision
// (i.e. an 80 or 128 bit double). Extended precision uses more memory ( 
// storing an extended precision double probably costs uses 16 bytes ), and more
// cpu time
//
// valid values for HIGHER_PRECISION include
//         "long double"
//         "__float128" ( for gcc )
//
// Note that the implementation of long double is compiler specific. 