.PHONY: default objects executable all clean
CC = gcc

LD = $(CC)
LDFLAGS = -lm

SOURCE_C = $(wildcard *.c)
OBJECTS_C = $(patsubst %.c, %_c.o, $(SOURCE_C))

DEPS_DIR = ../../src
DEPS_OBJ = $(DEPS_DIR)/master_c.o

CFLAGS = -std=gnu99 -I$(DEPS_DIR) -O0 -g
EXECUTABLE = ./test_eval.e

default: all

objects: $(OBJECTS_C)

executable: $(EXECUTABLE)

all: deps objects executable

deps:
	$(MAKE) -C $(DEPS_DIR) objects

%_c.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

%.e: $(OBJECTS_C)
	$(LD) $^ -o $@ $(DEPS_OBJ) -lm

clean:
	rm -rfv $(OBJECTS_C) $(EXECUTABLE)
