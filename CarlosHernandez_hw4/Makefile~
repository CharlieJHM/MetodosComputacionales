anim.gif:anim34.png
	./animate.sh
anim34.png:momo.txt
	python Plots_Temperatura.py
momo.txt:a.out	
	./a.out > momo.txt
a.out:DifusionTemperatura.c
	gcc DifusionTemperatura.c

clean:
	rm *.png *.txt
