#*****************************************************************************
#* Configuration
#*****************************************************************************

LIBNAME=airshowers
LIBRARY_DIR = ./lib

# External Libraries:
EXT_LIBS       = $(shell root-config --libs) -lMinuit
EXT_CFLAGS     = $(shell root-config --cflags) -fPIC
EXT_LDFLAGS    = root-config --ldflags -fPIC

# Local Directory Structure:
INCLUDE_DIRS   = ./include 
OBJECT_DIR     = ./tmp
LIBRARY_DIR    = ./lib
SRC_DIR        = ./src
TEST_DIR       = ./test
BIN_DIR        = ./run/bin

# Targets:
CPP_SOURCES    = $(wildcard $(SRC_DIR)/*.cpp) 
LIB_OBJECTS    = $(CPP_SOURCES:$(SRC_DIR)/%.cpp=$(OBJECT_DIR)/%.o) 
LIB_SHARED     = $(LIBRARY_DIR)/lib$(LIBNAME).so
EXE_SOURCES    = $(wildcard $(TEST_DIR)/*.cpp) 
EXECS          = $(EXE_SOURCES:$(TEST_DIR)/%.cpp=$(BIN_DIR)/%)

# Dependencies:
INCLUDE_DEPS   = $(wildcard ./include/*.h)

# Compiler Options:
CFLAGS         = $(INCLUDE_DIRS:%=-I%) $(EXT_CFLAGS) -O0
LDFLAGS        = $(LIBRARY_DIR:%=-L%) $(EXT_CFLAGS) -O0
CC             = g++
LD_STATIC      = ar cqs
LD_SHARED      = g++ -g --shared
RM             = rm -f

#*****************************************************************************
#* Dependencies
#*****************************************************************************
#all: install_dirs $(LIB_STATIC) $(LIB_SHARED) $(EXECS)

all: $(LIB_SHARED) $(EXECS)

# install_dirs:
#	@ install/install_dirs.csh

# $(LIB_STATIC): $(LIB_OBJECTS)
#	$(RM) $(LIB_STATIC)
#	$(LD_STATIC) $(LIB_STATIC) $(LIB_OBJECTS) $(EXT_LIBS)  

$(LIB_SHARED): $(LIB_OBJECTS)
	$(LD_SHARED) $(LIB_OBJECTS) -o $(LIB_SHARED) $(EXT_LIBS) $(EXT_SHARED)  

$(LIB_OBJECTS): $(INCLUDE_DEPS) Makefile

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(EXECS) : $(INCLUDE_DEPS)

$(BIN_DIR)/%: $(TEST_DIR)/%.cpp $(LIB_SHARED)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ $(EXT_LIBS) $(LIB_SHARED) $(EXT_SHARED)  

clean:
	$(RM) $(LIB_OBJECTS) $(EXT_OBJECTS)

binclean:
	$(RM) $(LIB_SHARED) $(EXECS)
	$(RM) $(LIB_STATIC) 
