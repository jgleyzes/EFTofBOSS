
theoryCovness=0
simname=LightConeHector
grid=LightConeHectorPatchyWideHDvFFT_IR06

cd /exports/pierre/EFTofBOSS

for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 1 2
do
qsub -l nodes=node01:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.Prior$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.Prior$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 
done
done
done
done
done
done

#################

for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 3 4
do
qsub -l nodes=node02:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 
done
done

done

done
done
done

#######################

for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 5 6
do
qsub -l nodes=node03:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done


#######################

for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 7 8
do
qsub -l nodes=node04:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done


#####################
#####################

for planck in 0
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 9 10
do
qsub -l nodes=node05:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done

#####################

for planck in 1
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 1 2
do
qsub -l nodes=node06:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done

#####################

for planck in 1
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 3 4
do
qsub -l nodes=node07:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done

#####################

for planck in 1
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 5 6
do
qsub -l nodes=node08:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done

#####################

for planck in 1
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 7 8
do
qsub -l nodes=node09:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done

#####################

for planck in 1
do
for prior in 20
do
for marg in 1
do
for kmax in 0.25
do
for kmaxbisp in 0.06 0.08 0.1
do
for number in 9 10
do
qsub -l nodes=node10:ppn=4 -N b$number.$kmax.$simname.$theoryCovness.$marg.$planck.$prior -o log/output_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -e log/error_b$number.$kmax.$kmaxbisp.$simname.Planck$planck.$prior.log -v aa=$number,b=$kmax,c=$simname,d=$planck,e=$grid,f=$kmaxbisp,g=$marg,h=$theoryCovness,ii=$prior cluster_MCMC_Mercury.sh 

done
done

done

done
done
done