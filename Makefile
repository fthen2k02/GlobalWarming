CXX = g++
CXXFLAGS = -Wall -std=c++17 -O3 -IAudioFile
TARGET = process_file.exe

all: $(TARGET)

$(TARGET): $*.cpp
	$(CXX) $(CXXFLAGS) -o $@ $**

clean:
	del $(TARGET)
