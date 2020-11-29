# Local folders
SRC ?= ./src
INC ?= ./include
OBJ ?= ./obj
BIN ?= ./bin

# Common binaries
GCC             ?= g++

#user flags
CCFLAGS   := -Wall -lfftw3 -lcfitsio

# Extra user flags
EXTRA_CCFLAGS = -I$(INC)

# Optimisation flags
CCFLAGS += -O3 

# DEBUG flags
#CCFLAGS   += -g

# PIC flags
#CCFLAGS   += -fPIC 

MAIN = $(SRC)/main.c
SOURCE_FILES = $(filter-out $(MAIN), $(wildcard $(SRC)/*.c))
OBJ_FILES = $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(SOURCE_FILES))

#####################
# default make rule #
#####################

TARGET = main

.PHONY: all
all: $(BIN)/$(TARGET)


$(BIN)/$(TARGET): $(MAIN) $(SOURCE_FILES)
	$(GCC) $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS)

$(OBJ)/%.o: $(SRC)/%.c
	$(GCC) -c $(CCFLAGS) $(EXTRA_CCFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -rf $(BIN)/$(TARGET) $(OBJ)/*
