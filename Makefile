complex?=complexes/coloring.txt # sonst make complex=complexes/complex.txt

run: build/Flag
	./build/Flag

simp: build/Simp
	./build/Simp $(complex)

build/Flag: build src/Flag.cc
	g++ -std=c++0x -o build/Flag src/Flag.cc

build/Simp: build src/Simp.cc
	g++ -std=c++0x -o build/Simp src/Simp.cc

build:
	mkdir build

clear:
	rm -rf build
