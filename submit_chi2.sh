
theoryCovness=0

simname=ChallengeA
grid=ChallengeHDvFFT_IR06

#simname=LightConeHector
#grid=LightConeHectorPatchyWideHDvFFT_IR06

#grid=LightConeHectorDataHDvFFT_IR06

cd /exports/pierre/EFTofBOSS



for planck in 0
do
for prior in 100
#for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0
do
for number in A
#for number in F
do
qsub -N c$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_chi2b$number.$kmax.$kmaxbisp.$simname.Planck$planck.Prior$prior.log -e log/error_chi2b$number.$kmax.$kmaxbisp.$simname.Planck$planck.Prior$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior chi2_mercury.sh 
done
done
done
done
done
done