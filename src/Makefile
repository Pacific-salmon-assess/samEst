all:  Ricker_simple_autocorr.so Ricker_tv_all.so SR_HMM_all.so  Rickerkf.so



Ricker_simple_autocorr.so: Ricker_simple_autocorr.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple_autocorr.cpp')"


SR_HMM_all.so: SR_HMM_all.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_all.cpp')"
		
Rickerkf.so: Rickerkf.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf.cpp')"

Ricker_tv_all.so: Ricker_tv_all.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tv_all.cpp')"




clean:
	rm -rf *o
