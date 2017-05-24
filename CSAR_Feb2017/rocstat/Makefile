
CC=gcc
GET=co

CFLAGS=-Wall -march=native -O3 -fopenmp
LDFLAGS=-lX11 -lXpm -lm

OBJS = amira.o help.o random.o slice.o stat.o hotspot.o view.o
SRCS = $(OBJS:.o=.c)

all: rocstat

$(SRCS):
	$(CC) $(CFLAGS) -c $*.c

rocstat: $(OBJS) rocstat.c
	${CC} $(CFLAGS) -o $@ $(OBJS) $@.c $(LDFLAGS)

install: rocstat
	install -d $(HOME)/bin
	install -C rocstat rocplot $(HOME)/bin

install-local: xview rocstat
	install -C rocstat rocplot /usr/local/bin

uninstall:
	$(RM) $(HOME)/bin/{rocstat,rocplot}

uninstall-local:
	$(RM) /usr/local/bin/{rocstat,rocplot}

clean:
	$(RM) $(OBJS) rocstat
