Makefile: 

PROJECT_DIR = /d/Git/N-Body-Simulation/experimental
OBJ_DIR = $(PROJECT_DIR)/objs

CC = g++
#LIBS = 

#INCLUDES =
#LIB_INCLUDES = 

#INC = $(INCLUDES) $(LIBINCLUDES)  

OBJS = main.o

main.o: main.cpp 
	$(CC) -c main.cpp
	mv main.o $(OBJ_DIR)/main.o


exe: clean $(OBJS)
	$(CC) -o $(PROJECT_DIR)/a.exe  $(OBJ_DIR)/*.o
 
clean: 
	rm -f $(OBJ_DIR)/*.o