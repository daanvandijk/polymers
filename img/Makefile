all: img

img: experiments-cm.png 

experiments-cm.png: experiments-cm.py
	./experiments-cm.py

clear: clean
clean:
	\rm *.png || true
	rm _pypolymer.so pypolymer.py || true
	rm __pycache__ -rf || true
