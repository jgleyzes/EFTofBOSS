
planck=0
theoryCovness=0
simname=LightConeDida
grid=LightConeDidaHDvFFT_IR06
#grid=LightConeHectorDataHDvCUBA_IR06


cd /exports/pierre/EFTofBOSS

#####################

#####################

for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.2 0.25 0.3
do
for kmaxbisp in 0
do
for number in data
do
qsub -l nodes=node09:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done

for planck in 1
do
for prior in 20
do
for marg in 1
do
for kmax in 0.2 0.25 0.3
do
for kmaxbisp in 0
do
for number in data
do
qsub -l nodes=node10:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done
