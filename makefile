#Compiler 
CC = g++

#include the header directory
INC_DIR = headers

#source file directory
SRC_DIR = src

#an object directory
OBJ_DIR = object

#flags for the compiler
CFLAGS = -std=c++11 -Wall -I.

SRCS = $(SRC_DIR)/heat1d.cpp  $(SRC_DIR)/heat2d.cpp  $(SRC_DIR)/testing.cpp  $(SRC_DIR)/vector_operations.cpp

OBJS = $(SRC_DIR)/heat1d.o  $(SRC_DIR)/heat2d.o  $(SRC_DIR)/testing.o  $(SRC_DIR)/vector_operations.o

DEPS = $(INC_DIR)/matrix.h $(INC_DIR)/vector.h $(INC_DIR)/heat1d.h $(INC_DIR)/heat2d.h   

all : $(OBJS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

$(OBJ_DIR)/conjugategradient.o : $(DEPS)

$(OBJ_DIR)/testing.o : $(DEPS)

clean:
	rm rf object
