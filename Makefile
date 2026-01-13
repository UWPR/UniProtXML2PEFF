all: UniProtXML2PEFF

UniProtXML2PEFF:
	g++ -O3 -std=c++11 -o UniProtXML2PEFF.exe UniProtXML2PEFF.cpp -ltinyxml2

clean:
	rm -f UniProtXML2PEFF.exe

test: all
	./UniProtXML2PEFF.exe test.xml test.peff
