#ifndef OOMPH_MATHEMATICA_DEFINITIONS
#define OOMPH_MATHEMATICA_DEFINITIONS

/* C language definitions for use with Mathematica output */
#define Power(x, y) (pow((double)(x), (double)(y)))
#define Sqrt(x)     (sqrt((double)(x)))

#define Abs(x)      (fabs((double)(x)))

#define Exp(x)      (exp((double)(x)))
#define Log(x)      (log((double)(x)))
#define Sin(x)      (sin((double)(x)))
#define Cos(x)      (cos((double)(x)))
#define Tan(x)      (tan((double)(x)))

#define ArcSin(x)       (asin((double)(x)))
#define ArcCos(x)       (acos((double)(x)))
//#define ArcTan(x)       (atan((double)(x)))
#define ArcTan(x,y)     (atan2((double)(y),(double)(x)))
#define ArcCot(x)       (atan(1.0/(double)(x)))
#define Sinh(x)         (sinh((double)(x)))
#define Cosh(x)         (cosh((double)(x)))
#define Tanh(x)         (tanh((double)(x)))


//#define E       2.71828182845904523536029
//#define Pi      3.14159265358979323846264
//#define Degree      0.01745329251994329576924

#endif
