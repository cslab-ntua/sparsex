#ifndef MACROS_H__
#define MACROS_H__

#define _CON3(a,b,c) a ## b ## c
#define CON3(a,b,c) _CON3(a,b,c)

#define _CON5(a,b,c,d,e) a ## b ## c ## d ## e
#define CON5(a,b,c,d,e) _CON5(a,b,c,d,e)

#define _CON6(a,b,c,d,e,f) a ## b ## c ## d ## e ## f
#define CON6(a,b,c,d,e,f) _CON6(a,b,c,d,e,f)

#define _CON7(a,b,c,d,e,f,g) a ## b ## c ## d ## e ## f ## g
#define CON7(a,b,c,d,e,f,g) _CON7(a,b,c,d,e,f,g)

#define _CON8(a,b,c,d,e,f,g,h) a ## b ## c ## d ## e ## f ## g ## h
#define CON8(a,b,c,d,e,f,g,h) _CON8(a,b,c,d,e,f,g,h)

#define UINT_TYPE(bits) CON3(uint, bits, _t)

#endif /* MACROS_H__ */
