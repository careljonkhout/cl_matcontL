#define SUNDIALS_DOUBLE_PRECISION 1
//#define SUNDIALS_EXTENDED_PRECISION 1
// comment line 1 and uncomment line 2 to enable extended precision 
// (i.e. an 80 bit double). Extended precision uses more memory ( storing a 
// extended precision double probably costs uses 16 bytes ), more cpu time
// note that extended precision, is defined in the C code as long double. The 
// implementation of long double is compiler specific. 