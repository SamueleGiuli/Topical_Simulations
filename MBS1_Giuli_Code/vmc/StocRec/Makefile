GFORT := gfortran

Task: ./vmc.x ./vmc_sr_read.x

vmc.x: ./VMC.f90
	$(GFORT) VMC.f90 -o vmc.x


vmc_sr_read.x: ./VMC_SR_read.f90
	$(GFORT) VMC_SR_read.f90 -o vmc_sr_read.x

clean:
	@rm *.x
	@rm modello.mod
