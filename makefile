#Compiler 
CC = g++

#include the header directory
INC_DIR = headers

#source file directory
SRC_DIR = src

#an object directory
OBJ_DIR = object

#directory for executables
EXE_DIR = executables

#directory that saves the data in txt files
DATA = plots

#flags for the compiler
CFLAGS = -std=c++11 -Iheaders

#SRCS = $(SRC_DIR)/heat1d.cpp  $(SRC_DIR)/heat2d.cpp  $(SRC_DIR)/testing.cpp 

OBJS = $(OBJ_DIR)/heat1d.o $(OBJ_DIR)/heat2d.o $(OBJ_DIR)/testing.o

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# The default rule
$(EXE_DIR)/testing: $(OBJS) | $(EXE_DIR)
	$(CC) $(CLFAGS) -o $@ $^

$(EXE_DIR):
	mkdir $@
	
$(OBJ_DIR):
	mkdir $@

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
	-rm -rf $(EXE_DIR)
	-rm $(DATA)/*txt
	-rm $(DATA)/*png
