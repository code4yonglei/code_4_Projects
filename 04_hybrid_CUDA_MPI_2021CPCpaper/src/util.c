#include <sys/time.h>
#include <sys/resource.h>

#include "util.h"

// The CPU time (in seconds)
// calls getrusage, limited accuracy

double second()
{
    static struct rusage temp;

    getrusage(RUSAGE_SELF,&temp);

    double foo1 = temp.ru_utime.tv_sec;  // seconds
    double foo2 = temp.ru_utime.tv_usec; // uSecs
    return  foo1 + (foo2/1000000.0); // milliseconds
}
