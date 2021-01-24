complex?=complexes/example.txt # sonst make complex=complexes/complex.txt

run: build/Simp
	./build/Simp $(complex)

build/Simp: build src/Simp.cc
	g++ -std=c++0x -o build/Simp src/Simp.cc

build:
	mkdir build

clear:
	rm -rf build
