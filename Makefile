CC = gcc
OBJECT = functions.o
EXECUTABLE = main

$(EXECUTABLE): $(OBJECT)
	$(CC) main.c $(OBJECT) -o $@ 

functions.o: functions.c
	$(CC) -c functions.c -o functions.o

clean:
	rm -f $(OBJECT) $(EXECUTABLE)
