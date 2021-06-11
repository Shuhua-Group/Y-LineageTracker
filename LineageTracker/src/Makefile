PRGS =  baseml codeml basemlg mcmctree pamp evolver yn00 chi2 
CC = cc # cc, gcc, cl

CFLAGS = -O3
#CFLAGS = -fast

LIBS = -lm # -lM

mcmctree : mcmctree.c  tools.c treesub.c treespace.c paml.h 
	$(CC) $(CFLAGS) -o $@ mcmctree.c tools.c $(LIBS)
