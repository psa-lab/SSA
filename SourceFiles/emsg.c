#include <stdio.h>
 
emsg(msg, source, sever, file, line)
char *msg;
int source, sever;
char *file;
int line;
{
#ifdef CMSWATC
        fflush(stderr);
        fflush(stdout);
        printf("%s %d", msg, source);
        fflush(stdout);
        sleep(2);
#else
        printf("In file %s, line %d: %s %d\n", file, line, msg, source);
#endif CMSWATC
}
