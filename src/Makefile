all:   Ricker_simple_autocorr.so  SR_HMM.so SR_HMM_a.so SR_HMM_b.so  Ricker_tv_all.so SR_HMM_all.so SR_HMM_all_logb.so



Ricker_simple_autocorr.so: Ricker_simple_autocorr.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple_autocorr.cpp')"

Ricker_autocorr_deprecated.so: Ricker_autocorr_deprecated.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_autocorr_deprecated.cpp')"

Ricker_simple_deprecated.so: Ricker_simple_deprecated.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_simple_deprecated.cpp')"

SR_HMM_all_logb.so: SR_HMM_all_logb.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_all_logb.cpp')"

SR_HMM_all.so: SR_HMM_all.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_all.cpp')"
		
SR_HMM.so: SR_HMM.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM.cpp')"

SR_HMM_a.so: SR_HMM_a.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_a.cpp')"

SR_HMM_b.so: SR_HMM_b.cpp
	Rscript --vanilla -e "TMB::compile('SR_HMM_b.cpp')"

Rickerkf.so: Rickerkf.cpp
	Rscript --vanilla -e "TMB::compile('Rickerkf.cpp')"

Ricker_tv_all.so: Ricker_tv_all.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tv_all.cpp')"

Ricker_tv_all_deprecated.so: Ricker_tv_all_deprecated.cpp
	Rscript --vanilla -e "TMB::compile('Ricker_tv_all_deprecated.cpp')"


clean:
	rm -rf *o
