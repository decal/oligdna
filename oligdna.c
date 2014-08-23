/*
   Generate Organic Quadruplet Nucleotide DNA Sequences From Base-N
   Oligodeoxyribonucleotide Primers for Polymerase Chain Reaction

   Written by Derek Callaway ( decal [AT] sdf [D0T] org )

   Compile: gcc -Wall -ansi -pedantic -o oligdna oligdna.c
*/

#include <stdio.h>
#include <stdlib.h>

void usage(const char *av0) {
  if(!av0)
    av0 = "specdna";

  fprintf(stderr, "usage: %s <file>\n", av0);

  exit(EXIT_FAILURE);
}

void vexit(const char *f) {
  if(!f)
    f = "vexit";

  perror(f);

  exit(EXIT_SUCCESS);
}

static char *_B="Tcg", *_D="ATg", *_H="ATc", *_K="Tg", *_M="Ac";
static char *_N="AcTg", *_R="Ag", *_S="cg", *_V="Acg", *_W="AT", *_Y="cT";

int main(int argc, char *argv[]) {
  char buf[BUFSIZ] = { 0 }, **mixs = NULL, **vals = NULL;
  register FILE *fp = NULL;
  register unsigned int acc = 0, cnt = 0;

  if(argc < 2)
    usage(*argv);

  fp = fopen(argv[1], "r");

  if(!fp)
    vexit("fopen");

  while(fgets(buf, sizeof buf, fp)) {
    register char *p = NULL, **pp = NULL, **pp2 = NULL;

    mixs = vals = NULL;

    for(p = buf, cnt = 0, acc = 0;*p;p++)
      switch(*p) {
        case 'B':
        case 'D':
        case 'H':
        case 'K':
        case 'M':
        case 'N':
        case 'R':
        case 'S':
        case 'V':
        case 'W':
        case 'Y':
          cnt++;
      }

      if(cnt >= acc) {
        register size_t alen = 1 + cnt;

        mixs = malloc(alen * sizeof *mixs);

        if(!mixs)
          vexit("malloc");

        vals = malloc(alen * sizeof *vals);

        if(!vals)
          vexit("malloc");
      }

      for(pp = mixs, pp2 = vals, p = buf;*p;p++)
        switch(*p) {
          case 'B':
            *pp++ = _B;
            *pp2++ = _B;

            break;
          case 'D':
            *pp++ = _D;
            *pp2++ = _D;

            break;
          case 'H':
            *pp++ = _H;
            *pp2++ = _H;

            break;
          case 'K':
            *pp++ = _K;
            *pp2++ = _K;

            break;
          case 'M':
            *pp++ = _M;
            *pp2++ = _M;

            break;
          case 'N':
            *pp++ = _N;
            *pp2++ = _N;

            break;
          case 'R':
            *pp++ = _R;
            *pp2++ = _R;

            break;
          case 'S':
            *pp++ = _S;
            *pp2++ = _S;

            break;
          case 'V':
            *pp++ = _V;
            *pp2++ = _V;

            break;
          case 'W':
            *pp++ = _W;
            *pp2++ = _W;

            break;
          case 'Y':
            *pp++ = _Y;
            *pp2++ = _Y;

            break;
          default:
            *pp++ = p;
            *pp2++ = p;
        }

        *pp = NULL;

        if(cnt)
          acc = --cnt;

        for(;*vals[acc];++*vals) {
          register unsigned int j = 0;

          if(!**vals) {
            do {
              vals[j] = mixs[j];
              ++vals[++j];
            } while(j < acc && !*vals[j]);

            if(j == acc && !*vals[acc])
              break;
           }

           for(j = 0, p = buf;*p;p++)
             switch(*p) {
               case 'B':
               case 'D':
               case 'H':
               case 'K':
               case 'M':
               case 'N':
               case 'R':
               case 'S':
               case 'V':
               case 'W':
               case 'Y':
                 putchar(*vals[j++]);

                 break;
               default:
                 putchar(*p);
             }
        }
  }

  exit(EXIT_SUCCESS);
}
