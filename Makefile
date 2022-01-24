CC = g++
TARGET = HTLP

SRC = $(wildcard *.cpp src/*.cpp)
FILE = $(notdir $(SRC))
OBJ = $(patsubst %.cpp, build/%.o, $(FILE))

# all:
# 	@echo $(SRC)
# 	@echo $(FILE)
# 	@echo $(OBJ)
# 	@echo "end"

FLAGES = -g -O2 -std=c++2a -pthread -march=native -c
LIB = -lntl -lssl -lcrypto -lgmp -lm

$(TARGET) : $(OBJ)
	$(CC) $^ -o $(TARGET) $(LIB)

build/%.o: src/%.cpp
	$(CC) $(FLAGES) $< -o $@

.PHONY:clean
clean:
	rm -f *.o src/*.o build/*.o $(TARGET)


	