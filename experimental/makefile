Makefile: 

#PROJECT_DIR = /cygdrive/d/Git/N-Body-Simulation/experimental
#OBJ_DIR = $(PROJECT_DIR)/objs

CC = g++
CC_FLAGS = -Dgtest_disable_pthreads=ON
#LIBS = 

#INCLUDES =
#LIB_INCLUDES = 

#INC = $(INCLUDES) $(LIBINCLUDES)  

OBJS = main.o

main.o: main.cpp 
	$(CC) $(CC_FLAGS) -c main.cpp
	mv main.o objs/main.o
	

main: clean # $(OBJS)
#	$(CC) $(CC_FLAGS) -o main.exe  objs/*.o
	$(CC) $(CC_FLAGS) -o main.exe main.cpp
clean: 
	rm -f $(OBJ_DIR)/*.o