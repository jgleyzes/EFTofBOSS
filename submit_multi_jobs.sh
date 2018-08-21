#submit z=0.61 jobs first, simtype is LightConeDida
#Simnames 

#simname=ChallengeQuarter_A
simname=LightConeHector

#for kmax in 0.1
#for kmax in 0.2 0.25
for kmax in 0.25
do
#for kmaxbisp in 0.07 0.08 
for kmaxbisp in 0.07
do
#for number in A B G F
for number in {1..9}
do
 #the above is for challenge boxes
 #sbatch --job-name=b$number.$kmax.$simname --output=outfiles/b$number.$kmax.$simname.out cluster_MCMC.sh $number $kmax $simname\_$number 0 conLightConeHectorwideh
sbatch --job-name=b$number.$kmax.$simname --output=outfiles/b$number.$kmax.$kmaxbisp.$simname.out cluster_MCMC.sh $number $kmax $simname 0 conLightConeHectorwideh $kmaxbisp
done
sbatch --job-name=bdata.$kmax.$simname --output=outfiles/bdata.$kmax.$kmaxbisp.$simname.out cluster_MCMC.sh data $kmax $simname 0 conLightConeHectorwideh $kmaxbisp
done
done

