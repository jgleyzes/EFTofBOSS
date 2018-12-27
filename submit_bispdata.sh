theoryCovness=0
simname2=LightConeHector
grid=LightConeHectorDataHDvFFT_IR06

#################

cd /exports/pierre/EFTofBOSS


for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25 0.3
do
for kmaxbisp in 0.08 0.1
do
for number in data
do
qsub -l nodes=node07:ppn=4 -N b$number.$kmax.$simname2.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname2.Planck$planck.Prior$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname2.Planck$planck.Prior$prior.log -v aa=$number,b=$kmax,c=$simname2,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 
done
done
done
done
done
done

#################

for planck in 1
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25 0.3
do
for kmaxbisp in 0.08 0.1
do
for number in data
do
qsub -l nodes=node06:ppn=4 -N b$number.$kmax.$simname2.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname2.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname2.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname2,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 
done
done

done

done
done
done