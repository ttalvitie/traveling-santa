.PHONY: clean

all: optimize

optimize: optimize.cpp
	g++ optimize.cpp -o optimize -O2 -Wall -Wextra -std=c++17 -g

clean:
	rm -f optimize
