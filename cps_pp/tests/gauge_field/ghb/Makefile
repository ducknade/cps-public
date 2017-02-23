all :

qcd :
	make clean
	make -C binaries
	mkdir results
	cd scripts && time ../binaries/NOARCH.x | tee ../results/log

bgq :
	make clean
	make -C binaries
	mkdir results
	cd scripts && run-job.sh ./run.sh | tee ../results/log

clean :
	make -C binaries clean || true
	rm -rf results || true
