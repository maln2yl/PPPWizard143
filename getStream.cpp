#include <cstdio>
#include <unistd.h>
#include <vector>
#include "rtklib.h"

typedef struct
{
  stream_t stream;
  int source;
  int format;
} STREAM;

int showmsg(char *format,...)
{
  va_list arg;
  va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
  fprintf(stderr,*format?"\r":"\n");
  return 0;
}

void settspan(gtime_t ts, gtime_t te)
{
}

void settime(gtime_t time)
{
}

int main(int argc, char *argv[])
{
  STREAM ls;
  std::vector<STREAM> stream;

  while (!feof(stdin))
  {
    char s[100], str[100];
    *s='\0';
    //gets(s);
    fgets(s, sizeof s, stdin);
    if (*s)
    {
//    printf("%s\n", s);
      strinit(&ls.stream);
      sscanf(s, "%s%d%d", str, &ls.source, &ls.format);
      stropen(&ls.stream, ls.source, STR_MODE_R, str);
      stream.push_back(ls);
    }
  }
  
  while (1)
  {
    int n;
    gtime_t time=timeget();
    for (n=0;n<(int)stream.size();n++)
    {
      int i;
      unsigned char buff[50];
      i=strread(&stream[n].stream, buff, 50);
      if (i)
      {
        int j;
        printf("%d %d %ld ", n+1, stream[n].format, time.time);
        for (j=0;j<i;j++)
        {
          printf("%02X", buff[j]);
        }
        printf("\n");
        fflush(stdout);
      }
    }
    usleep(10000);
  }

  return 1;
}
