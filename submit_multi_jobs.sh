#submit z=0.61 jobs first, simtype is LightConeDida
#Simnames 

simname=ChallengeA
#simname=LightConeHector

#for kmax in 0.1
for kmax in 0.2 0.25 0.3
#for kmax in 0.25
do
for kmaxbisp in 0.07 0.08 
#for kmaxbisp in 0.07
do
for number in A B G F
#for number in {1..9}
do
 #the above is for challenge boxes
sbatch --job-name=b$number.$kmax.$simname --output=outfiles/b$number.$kmax.$kmaxbisp.$simname.marg.out cluster_MCMC.sh $number $kmax $simname 0 conLightConeHectorwideh $kmaxbisp
done
#sbatch --job-name=bdata.$kmax.$simname --output=outfiles/bdata.$kmax.$kmaxbisp.$simname.marg.out cluster_MCMC.sh data $kmax $simname 0 conLightConeHectorwideh $kmaxbisp
done
done

